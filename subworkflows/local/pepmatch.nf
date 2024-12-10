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


process filterResults {
    container = 'jbrestel/iedb'

    input:
    path(csv)

    output:
    path("filtered.csv"), emit: csv

    script:
    """
    awk -F, 'NR > 1 && \$2 != "" { OFS=","; print \$1, \$3, \$8, \$10, \$11 }' $csv > filtered.csv
    """

}

workflow pepMatch {
    take:
    peptideFasta
    proteomeFasta

    main:
    preprocess(proteomeFasta)
    pepMatchResults = match(peptideFasta, preprocess.out)
    filteredResults = filterResults(pepMatchResults.csv)

    emit:
    filteredResults.csv
}
