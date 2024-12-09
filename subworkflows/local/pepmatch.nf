#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process preprocess {
//    container = 'veupathdb/epitopemapping'
    container = 'jbrestel/iedb'
    input:
    path(fasta)

    output:
    path("reference_proteome*")

    script:
    """
    pepmatch-preprocess -p $fasta -k ${task.ext.kmer_size} -f ${task.ext.format} -n reference_proteome
    """
}

process match {
    //container = 'veupathdb/epitopemapping'
    container = 'jbrestel/iedb'

    input:
    path(queryFasta)
    path(preprocessedResults)

    output:
    path("pepMatchResults.csv"), emit: csv

    script:
    """
    pepmatch-match -q $queryFasta -p reference_proteome -k ${task.ext.kmer_size} -m ${task.ext.num_mismatches} -o pepMatchResults.csv
    """
}


workflow pepMatch {
    take:
    peptideFasta
    proteomeFasta

    main:
    preprocess(proteomeFasta)
    pepMatchResults = match(peptideFasta, preprocess.out)

    emit:
    pepMatchResults.csv
}
