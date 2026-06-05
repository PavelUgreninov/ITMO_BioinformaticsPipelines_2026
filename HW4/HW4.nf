#!/usr/bin/env nextflow

include { BCFTOOLS_MPILEUP as MPILEUP } from './modules/nf-core/bcftools/mpileup/'
include { RUN_QC } from './modules/run_qc/'
include { TRIMM } from './modules/trimm/'
include { ASSEMBLE } from './modules/assemble/'
include { MAP_AND_PLOT } from './modules/map_and_plot/'
include { DOWNLOAD_SRA} from './modules/download_sra/'
include { VCF_FILTER } from './modules/vcf_filter/'

workflow trimmQc {

    take: 
        trimmed_data
    main:
        RUN_QC('trimmed', trimmed_data)
    emit:
        RUN_QC.out
}

workflow assessData {
    //This named workflow is developed to extract information either from SRA or local files
    take:
    input_stream

    main:
    input_stream.branch { meta, file_or_sra ->
        sra: file_or_sra =~ /^(SRR|ERR|DRR|SRX|ERX|DRX|SRS|ERS|DRS)[0-9]+$/
        local_path: true
    }
    .set { ch_split }

    ch_downloaded = DOWNLOAD_SRA(ch_split.sra)

    raw_sra = ch_downloaded
        .map { meta, files ->
            def new_meta = meta + [ single_end: files.size() == 1 ]
            return [ new_meta, files.sort() ]
        }

    raw_local = ch_split.local_path
        .map { meta, pathStr ->
            def files = []
            if (pathStr instanceof List) {
                files = pathStr.collect { file(it) }
            } else {
                files = [file(pathStr)]
            }
            return [meta, files]
        }
        .map { meta, files ->
            def new_meta = meta + [ single_end: files.size() == 1 ]
            return [ new_meta, files ]
        }

    raw_reads = raw_sra.mix(raw_local)

    emit:
    data = raw_reads
}

workflow {
    //If a csv file is provided
    if (params.csv_file) {
        ch_samplesheet = Channel.fromPath(params.csv_file)
            .splitCsv(header: true, strip: true, quote: '"')
            .map { row ->
                def meta = [id: row.id, type: row.type, repeat: row.repeat]
                def files = row.file_or_sra?.trim()?.split(',')?.collect { it.trim() }
                def file_or_sra = files.size() == 1 ? files[0] : files
                if (meta.type == 'reference') file_or_sra = file(file_or_sra)
                return [meta, file_or_sra]
            }
            .branch { meta, file_or_sra ->
                reference: meta.type == 'reference'
                sample:    meta.type == 'sample'
            }
            .set { ch_split }

        ch_references = ch_split.reference
        ch_samples    = ch_split.sample
        raw_reads = assessData(ch_samples).data
        ch_references_keyed = ch_references.map { meta, file -> [meta.id, file] }
    }
    //If a csv file is not provided
    else {
        if (params.numberSRA) {
            def sra_list = params.numberSRA.split(',')*.trim()
            ch_input = Channel.fromList(sra_list)
                .map { sra_id ->
                    def meta = [id: sra_id, type: 'sample', repeat: 1]
                    return [meta, sra_id]
                }
        } else if (params.input_reads_folder) {
            ch_input = Channel.fromPath("${params.input_reads_folder}/*.{fq,fastq,fastq.gz,fq.gz}")
                .map { file ->
                    def id = file.name.toString().tokenize('_').first()
                    return [id, file]
                }
                .groupTuple()
                .map { id, files ->
                    def meta = [id: id, type: 'sample', repeat: 1]
                    return [meta, files.sort()]
                }
        } else {
            error "Either --csv_file, --numberSRA, or --input_reads_folder must be provided"
        }

        raw_reads = assessData(ch_input).data

        //if allignment is to be performed
        if (params.ref_type == 'alignment') {
            if (!params.reference) error "For ref_type 'alignment' you must provide --reference FASTA"
            ch_references = Channel.of([ [id: 'ref'], file(params.reference) ])
        //if assembly is to be performed
        } else if (params.ref_type == 'assembly') {
            RUN_QC('raw', raw_reads)
            ch_trimmed_raw = TRIMM(raw_reads, params.adapters)
            ch_trimmed_qc = trimmQc(ch_trimmed_raw)
            ch_all_reads = ch_trimmed_raw
                .map { meta, files -> files }
                .collect()
                .map { listOfLists -> listOfLists.flatten() }
                .map { all_files -> file(all_files) }
            ASSEMBLE(ch_all_reads)
            ch_references = ASSEMBLE.out.reference.map { ref -> [ [id: 'ref'], ref ] }
        } else {
            error "params.ref_type must be either 'alignment' or 'assembly'"
        }
        ch_references_keyed = ch_references.map { meta, file -> [meta.id, file] }
    }

    if (!params.csv_file && params.ref_type == 'assembly') {
        ch_trimmed_for_mapping = ch_trimmed_raw
    } else {
        RUN_QC('raw', raw_reads)
        ch_trimmed_raw = TRIMM(raw_reads, params.adapters)
        ch_trimmed_qc = trimmQc(ch_trimmed_raw)
        ch_trimmed_for_mapping = ch_trimmed_raw
    }

    ch_trimmed_processed = ch_trimmed_for_mapping
            .map { meta, files ->
                def trimmed_reads = files.findAll { it.toString().contains("_p.fq.gz") }
                if (trimmed_reads.isEmpty()) trimmed_reads = files
                [meta, trimmed_reads.sort()]
            }
            .map { meta, reads ->
                def clean_meta = meta.subMap(['type', 'repeat'])
                [meta.id, clean_meta, reads]
            }

    ch_combined = ch_trimmed_processed.combine(ch_references_keyed, by: 0)
    ch_for_mapping = ch_combined.map { id, clean_meta, reads, ref_file ->
        def meta_with_id = clean_meta + [id: id]
        [meta_with_id, ref_file, reads]
    }

    MAP_AND_PLOT(ch_for_mapping)
    ch_mpileup_input = MAP_AND_PLOT.out.bam_out.map { meta, bam, bed1, bed2 -> [meta, bam, [], []] }
    ch_ref_input = MAP_AND_PLOT.out.ref_out.map { meta, ref, fai -> [meta, ref, fai] }
    MPILEUP(ch_mpileup_input, ch_ref_input, Channel.value(false))

    ch_filter_input = MPILEUP.out.vcf.join(MPILEUP.out.tbi).map { meta, vcf, tbi -> [meta, vcf, tbi] }
    VCF_FILTER(ch_filter_input)

    ch_final_vcf = VCF_FILTER.out.vcf
}