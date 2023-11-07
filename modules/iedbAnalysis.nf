#!/usr/bin/env nextflow
nextflow.enable.dsl=2 

process downloadIedbEpitope {

    input:
    path weblink

    output:

    script: 

    """
    wget ${weblink}
    """
}

workflow iedp {

    take:
    seq

    main:

    downloadSequences = downloadIedbEpitope(seq)
}