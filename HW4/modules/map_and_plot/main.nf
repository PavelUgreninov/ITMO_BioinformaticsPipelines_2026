process MAP_AND_PLOT {
    tag "$meta.id"

    publishDir "${params.outdir}/analysis", mode: 'copy'

    input:
        tuple val(meta), path(ref), path(reads)
    
    output:
        tuple val(meta), path("${meta.id}.sorted.bam"), path("${params.bed_file}"), path("${params.bed_file2}"), emit: bam_out
        tuple val(meta), path("${ref}"), path("${ref}.fai"), emit: ref_out

    script:
    def reads_cmd = reads instanceof List ? reads.join(' ') : reads
    """
    bwa index $ref
    bwa mem $ref $reads_cmd | samtools view -Sb - | samtools sort -o ${meta.id}.sorted.bam -
    samtools index ${meta.id}.sorted.bam
    samtools faidx $ref
    
    samtools depth ${meta.id}.sorted.bam > depth.txt
    
    python3 <<CODE
import matplotlib.pyplot as plt
import pandas as pd
import os

if os.path.exists('depth.txt') and os.path.getsize('depth.txt') > 0:
    d = pd.read_csv('depth.txt', sep='\\t', header=None)
    plt.figure(figsize=(10,6))
    plt.plot(d[2], color='blue')
    plt.fill_between(range(len(d[2])), d[2], color='blue', alpha=0.3)
    plt.title('Coverage Profile: ${meta.id}')
    plt.xlabel('Genomic Position')
    plt.ylabel('Read Depth')
    plt.savefig('${meta.id}_coverage.png')
else:
    plt.figure()
    plt.text(0.5, 0.5, 'No Coverage Data', ha='center')
    plt.savefig('${meta.id}_coverage.png')
CODE
    """
}