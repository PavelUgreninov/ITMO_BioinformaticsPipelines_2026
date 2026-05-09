# ITMO_BioinformaticsPipilines_2026

This pipelines allows you to call variants.
Standard usage: 
* nextflow run HW3.nf -profile container --input_reads_folder <fodler_with_reads> --ref_type alignment --reference <reference_fasta>
* nextflow run HW3.nf -profile container --input_reads_folder <fodler_with_reads> --ref_type aseemble
* using SRA number for retrieving files is also allowed
Profiles:
1. Test1: assemble
2. Test2: alignment
3. Local: uses conda environment on the local machine
4. Cluster: uses .yaml file with dependencies for conda environment
5. Container: uses docker image from DocHub to run processes

