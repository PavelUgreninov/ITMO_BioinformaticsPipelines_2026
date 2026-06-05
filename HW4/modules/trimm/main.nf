process TRIMM {
    tag "${meta.id}${meta.repeat ? '_rep' + meta.repeat : ''}"

    input:
    tuple val(meta), path(reads)
    path adapters

    output:
    tuple val(meta), path("out_${meta.id}*.fq.gz"), emit: trimmed_pairs

    script:
    def rep_suffix = meta.repeat ? "_rep${meta.repeat}" : ""
    def reads_list = reads instanceof List ? reads : [reads]
    
    if (reads_list.size() == 1) {
        """
        trimmomatic SE \\
            ${reads_list[0]} \\
            out_${meta.id}${rep_suffix}_p.fq.gz \\
            ILLUMINACLIP:${adapters}:2:30:10 LEADING:3 TRAILING:3 MINLEN:36
        """
    } else {
        """
        trimmomatic PE \\
            ${reads_list[0]} ${reads_list[1]} \\
            out_${meta.id}${rep_suffix}_R1_p.fq.gz out_${meta.id}${rep_suffix}_R1_u.fq.gz \\
            out_${meta.id}${rep_suffix}_R2_p.fq.gz out_${meta.id}${rep_suffix}_R2_u.fq.gz \\
            ILLUMINACLIP:${adapters}:2:30:10 LEADING:3 TRAILING:3 MINLEN:36
        """
    }
}