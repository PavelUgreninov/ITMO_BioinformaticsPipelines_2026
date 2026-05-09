#!/usr/bin/env nextflow

include { BCFTOOLS_MPILEUP as mpileup } from './modules/nf-core/bcftools/mpileup/'
include { run_qc } from './modules/run_qc/'
include { trimm } from './modules/trimm/'
include { assemble } from './modules/assemble/'
include { map_and_plot } from './modules/map_and_plot/'

workflow trimm_qc {

    take: 
        trimmed_data
    main:
        run_qc('trimmed', trimmed_data)
    emit:
        run_qc.out
}

workflow {
    if (params.numberSRA) {
        raw_reads = Channel.fromSRA(params.numberSRA)
            .map { id, files ->
                    def meta = [:]
                    meta.id = id
                    meta.single_end = files.size() == 1 ? true : false
                    return [ meta, files.sort() ]
                }
    } else {
        raw_reads = Channel.fromPath("${params.input_reads_folder}/*.{fq,fastq,fastq.gz,fq.gz}")
            .map { file ->
                def id = file.name.toString().tokenize('_').first().replaceFirst(/\.fastq(\.gz)?$/, '').replaceFirst(/\.fq(\.gz)?$/, '')
                return [ id, file ]
            }
            .groupTuple()
            .map { id, files ->
                def meta = [:]
                meta.id = id
                meta.single_end = files.size() == 1 ? true : false
                return [ meta, files.sort() ]
            }
    }

    run_qc('raw', raw_reads)
    trimmed_reads = trimm(raw_reads, params.adapters)
    trimm_qc(trimmed_reads.collect())

    if (params.ref_type == 'alignment') {
        reference_genome = file(params.reference, checkIfExists: true)
    } else {
        reference_genome = assemble(trimmed_reads.collect())
    }
    
    alignment_results = map_and_plot(reference_genome, trimmed_reads.collect())

    mpileup(map_and_plot.out.bam_out, map_and_plot.out.ref_out, Channel.value(false))

}