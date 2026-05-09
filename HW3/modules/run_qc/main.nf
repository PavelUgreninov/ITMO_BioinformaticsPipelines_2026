process run_qc {
    
    tag "$meta.id"

    publishDir { "${params.outdir}/${reads_type}_qc" }, mode: 'copy'

    input:
        val reads_type
        tuple val(meta), path(reads)
    output:
        path "${reads_type}_qc_report/", type: 'dir'

    script:
    """
    mkdir -p ${reads_type}_qc_report
    fastqc -o ${reads_type}_qc_report/ $reads
    """
}