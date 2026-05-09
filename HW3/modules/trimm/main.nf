process trimm {
    
    tag "${meta.id}"

    input:
        tuple val(meta), path(reads)
        path adapters

    output:
        tuple val(meta), path("out_${meta.id}_R*_p.fq.gz"), emit: trimmed_pairs

    script:
    if (meta.single_end) {
    """
    trimmomatic SE \\
            ${reads} \\
            out_${meta.id}_p.fq.gz \\
            ILLUMINACLIP:${adapters}:2:30:10 LEADING:3 TRAILING:3 MINLEN:36
    """
    } else {
    """
    trimmomatic PE ${reads[0]} ${reads[1]} \\
        out_${meta.id}_R1_p.fq.gz out_${meta.id}_R1_u.fq.gz \\
        out_${meta.id}_R2_p.fq.gz out_${meta.id}_R2_u.fq.gz \\
        ILLUMINACLIP:${adapters}:2:30:10 LEADING:3 TRAILING:3 MINLEN:36
    """
    }
}