#!/usr/bin/env nextflow
nextflow.enable.dsl=2 

process downloadSequences{
    
    input:
     path(url)

    output:
     path('iedb_export/*'), emit: xmlFielPath

    script:
     template 'downloadIedb.bash'


}


process processIedbEpitope {

    input:
     path(xmlFile) 
 
    output:
     path("*fa")

    script: 
     sample_base = xmlFile.getSimpleName()
     template 'generateIedbPeptidefasta.bash' 
}

workflow iedp {

    take:
    url

    main:
    database = downloadSequences(url)
    fileList = database.xmlFielPath.view()
    processSequences = processIedbEpitope(database.fileList)
}