#!/usr/bin/env nextflow

params.input              = '' 
params.input_reads_folder = './data/'
params.is_sra             = false
params.ref_type           = 'assemble' 
params.reference          = 'data/ref.fa'
params.outdir             = 'results'
params.adapters           = "${projectDir}/adapters/TruSeq3-PE.fa"

process run_qc {
    tag "$reads_label"
    publishDir "${params.outdir}/${reads_type}_qc", mode: 'copy'
    conda 'bioconda::fastqc=0.11.9'

    input:
        val reads_type
        tuple val(reads_label), path(reads)
    output:
        path "${reads_type}_qc_report/", type: 'dir'

    script:
    """
    mkdir ${reads_type}_qc_report
    fastqc -o ${reads_type}_qc_report/ $reads
    """
}

process trimm {
    tag "$reads_label"
    conda 'bioconda::trimmomatic=0.39' 
    
    input:
        tuple val(reads_label), path(reads)
    output:
        tuple val(reads_label), path("out_${reads_label}_R*_p.fq.gz")

    script:
    """
    trimmomatic PE ${reads[0]} ${reads[1]} \\
        out_${reads_label}_R1_p.fq.gz out_${reads_label}_R1_u.fq.gz \\
        out_${reads_label}_R2_p.fq.gz out_${reads_label}_R2_u.fq.gz \\
        ILLUMINACLIP:${params.adapters}:2:30:10 LEADING:3 TRAILING:3 MINLEN:36
    """
}

process assemble {
    tag "$reads_label"
    conda 'bioconda::spades=3.15.5 python=3.11'
    
    input: 
        tuple val(reads_label), path(reads)
    output: 
        path("scaffolds.fasta")

    script:
    """
    spades.py -1 ${reads[0]} -2 ${reads[1]} -o spades_out
    mv spades_out/scaffolds.fasta .
    """
}

process map_and_plot {
    tag "$reads_label"
    conda 'bioconda::bwa=0.7.17 bioconda::samtools=1.15 conda-forge::matplotlib=3.5 conda-forge::pandas=1.4'
    publishDir "${params.outdir}/analysis", mode: 'copy'

    input:
        path ref
        tuple val(reads_label), path(reads)
    
    output:
        path "${reads_label}_coverage.png"
        path "${reads_label}.sorted.bam"

    script:
    """
    bwa index $ref
    bwa mem $ref $reads | samtools view -Sb - | samtools sort -o ${reads_label}.sorted.bam -
    samtools index ${reads_label}.sorted.bam
    
    samtools depth ${reads_label}.sorted.bam > depth.txt
    
    python3 -c "
    import matplotlib.pyplot as plt
    import pandas as pd
    import os

    if os.path.exists('depth.txt') and os.path.getsize('depth.txt') > 0:
        d = pd.read_csv('depth.txt', sep='\\t', header=None)
        plt.figure(figsize=(10,6))
        plt.plot(d[2], color='blue')
        plt.fill_between(range(len(d[2])), d[2], color='blue', alpha=0.3)
        plt.title('Coverage Profile: ${reads_label}')
        plt.xlabel('Genomic Position')
        plt.ylabel('Read Depth')
        plt.savefig('${reads_label}_coverage.png')
    else:
        plt.figure()
        plt.text(0.5, 0.5, 'No Coverage Data', ha='center')
        plt.savefig('${reads_label}_coverage.png')
    "
    """
}

workflow trimm_qc {
    take: trimmed_data
    main:
        run_qc('trimmed', trimmed_data)
    emit:
        run_qc.out
}

workflow {
    if (params.is_sra) {
        raw_reads = Channel.fromSRA(params.input)
            .map { meta, files -> [ meta.id, files ] }
    } else {
        raw_reads = Channel.fromFilePairs("${params.input_reads_folder}/*_{1,2}.fq", checkIfExists: true)
    }

    run_qc('initial', raw_reads)
    trimmed_reads = trimm(raw_reads)
    trimm_qc(trimmed_reads)

    // Reference logic
    if (params.ref_type == 'input') {
        reference_genome = file(params.reference, checkIfExists: true)
    } else {
        // Assembles the first sample to use as a reference for all mapping
        reference_genome = assemble(trimmed_reads.first())
    }

    // Final Analysis
    map_and_plot(reference_genome, trimmed_reads)
}