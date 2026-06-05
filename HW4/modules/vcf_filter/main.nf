process VCF_FILTER {
    tag "${meta.id}_rep${meta.repeat}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'}"

    input:
        tuple val(meta), path(vcf), path(tbi)

    output:
        tuple val(meta), path("${meta.id}_rep${meta.repeat}.filtered.vcf.gz"), emit: vcf
        tuple val(meta), path("${meta.id}_rep${meta.repeat}.filtered.vcf.gz.tbi"), emit: tbi

    script:
    def prefix = "${meta.id}_rep${meta.repeat}"
    """
    bcftools filter \\
        -i 'QUAL > 20' \\
        -o ${prefix}.filtered.vcf.gz \\
        -Oz \\
        ${vcf}
    tabix ${prefix}.filtered.vcf.gz
    """

    stub:
    def prefix = "${meta.id}_rep${meta.repeat}"
    """
    touch ${prefix}.filtered.vcf.gz
    touch ${prefix}.filtered.vcf.gz.tbi
    """
}