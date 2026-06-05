process DOWNLOAD_SRA {
    tag "$sra_id"
    
    input:
    tuple val(meta), val(sra_id)
    
    output:
    tuple val(meta), path("*.fastq.gz")
    
    script:
    """
    fasterq-dump ${sra_id} --split-files
    gzip *.fastq
    """
}