process ASSEMBLE {

    tag "$meta.id"
    
    input: 
        tuple val(meta), path(reads)
    output: 
        path("scaffolds.fasta")

    script:
    """
    spades.py -1 ${reads[0]} -2 ${reads[1]} -o spades_out
    mv spades_out/scaffolds.fasta .
    """
}
