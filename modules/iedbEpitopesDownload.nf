#!/usr/bin/env nextflow
nextflow.enable.dsl=2 

process downloadSequences{
    
    publishDir "${params.results}/xmlFiles", mode: 'copy'

    input:
     path(url)

    output:
     path('test/*'), emit: xmlFielPath

    script:
     template 'downloadIedb.bash'


}


process processIedbEpitope {

    publishDir "${params.results}/IndividulaFasta", mode: 'copy'
    input:
     path(xmlFile) 
 
    output:
     path("*fa"), emit: fastaFile

    script: 
     sample_base = xmlFile.getSimpleName()
     template 'generateIedbPeptidefasta.bash' 
}


process mergeEpitopesFasta{

    publishDir "${params.results}/iedbFasta", mode: 'copy'
    input:
    path(fasta)

    output:
    path("*fa")
    script:

    """
    cat ${fasta} > iedbEpitpes.fa
    """
}


workflow iedbEpitopesDownload {

    take:
    url

    main:
    database = downloadSequences(url)
    fileList = database.xmlFielPath.flatten()
    processSequences = processIedbEpitope(fileList)
    allFasta = processSequences.fastaFile.collect()
    mergeFasta = mergeEpitopesFasta(allFasta)
}