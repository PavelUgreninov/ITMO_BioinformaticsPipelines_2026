# ITMO_BioinformaticsPipilines_2026

This pipelines allows you to call variants.
Standard usage:
* `nextflow run HW3.nf -profile container --input_reads_folder <fodler_with_reads> --ref_type alignment --reference <reference_fasta>`
* `nextflow run HW3.nf -profile container --input_reads_folder <fodler_with_reads> --ref_type aseemble`
* `using SRA number for retrieving files is also allowed`
Profiles:
1. Test1: test assemble
2. Test2: test alignment
3. Local: uses conda environment on the local machine
4. Cluster: uses .yaml file with dependencies for conda environment
5. Container: uses docker image from DocHub to run processes

To run this pipeline with summary in a CSV file:
* nextflow run HW4.nf --csv_file ./test_data/samplesheet_sra.csv -profile cluster
  or
  in the ./test_data/reads_hw4/ directory run:
```bash
fasterq-dump ERR16834488 --split-files
fasterq-dump ERR16834491 --split-files
fasterq-dump ERR16834499 --split-files
fasterq-dump SRR30011709 --split-files
fasterq-dump SRR30011711 --split-files
fasterq-dump SRR30011716 --split-files
```
and then run:
* nextflow run HW4.nf --csv_file ./test_data/samplesheet_files.csv -profile cluster
