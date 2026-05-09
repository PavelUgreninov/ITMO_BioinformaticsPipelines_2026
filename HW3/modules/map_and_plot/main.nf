process map_and_plot {
    tag "$meta.id"

    publishDir "${params.outdir}/analysis", mode: 'copy'

    input:
        path ref
        tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), path("${params.bed_file}"), path("${params.bed_file2}"), emit: bam_out
    tuple val(meta), path("${ref}"), path("${ref}.fai"), emit: ref_out


    script:
    """
    bwa index $ref
    bwa mem $ref $reads | samtools view -Sb - | samtools sort -o ${meta.id}.sorted.bam -
    samtools index ${meta.id}.sorted.bam
    samtools faidx $ref
    
    samtools depth ${meta.id}.sorted.bam > depth.txt
    
    python3 -c "
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
    "
    """
}
