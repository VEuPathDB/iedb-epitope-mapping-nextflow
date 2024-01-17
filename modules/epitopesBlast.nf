#!/usr/bin/env nextflow
nextflow.enable.dsl=2 


process peptideSimilarity {

    //container = 'veupathdb/epitopemapping'
    container = 'epitopemapping'

    publishDir "${params.results}", mode: 'copy'


    input:
      path(refFasta)
      path(pepfasta)
      path(pepTab)
      val(taxa)

    output:
      path("*Peptides.fasta"), emit: peptideFasta
      path("*PeptideGene.txt"), emit: pepResults

    script:
      template 'ProcessPeptides.bash'

}

process makeBlastDatabase {

     //container = 'veupathdb/blastsimilarity'
     container = 'epitopemapping'

    input:
      path(fasta)

    output:
      path("BlastDB"), emit: db
    script:
       sample_base = fasta.getSimpleName()
       template 'makeBlastDb.bash'

}


process blastSeq {

    //container = 'veupathdb/blastsimilarity'
    container = 'epitopemapping'

    publishDir "${params.results}/BlastOut", mode: 'copy'

    input:
      path(query)
      path(db)

    output:
      path("${sample_base}*xml"), emit: result

    script:
      sample_base = query.getSimpleName()
      template 'BlastSeq.bash'

}

process diamondDatabase {

    container = 'veupathdb/diamond'

    input:
      path(fasta)

    output:
      path("*dmnd"), emit: db

    script:
      sample_base = fasta.getSimpleName()
      template 'makeDiamondDb.bash'
    
}

process diamondBlast {

    container = 'veupathdb/diamond'

    publishDir "${params.results}/BlastOut", mode: 'copy'

    input:
      path(query)
      path(db)

    output:
      path("*xml"), emit: result

    script:
      sample_base = query.getSimpleName()
      template 'diamondBlast.bash'

}   

process processXml {

     //container = 'veupathdb/epitopemapping'
     container = 'epitopemapping'

    publishDir "${params.results}/BlastOut", mode: 'copy'

    input:
      path(xml)

    output:
      path("*txt"), emit: resultFormated


    script:
      template 'processBlastXml.bash'
}


process mergeeResultsFiles {

    //container = 'veupathdb/epitopemapping'
    container = 'epitopemapping'

    publishDir "${params.results}/BlastOut", mode: 'copy'

    input:
      path(exactMatch)
      path(balst)
    
  

    output:
      path("*txt")

    script:
      template 'mergeFiles.bash'

}

workflow epitopesBlast {

   take: 
    refFasta
    peptidesTab
    peptidesGeneFasta

    main:

    processPeptides = peptideSimilarity(refFasta, peptidesGeneFasta, peptidesTab, params.taxon)

    if (params.blastMethod == "ncbiBlast") {
      database = makeBlastDatabase(processPeptides.peptideFasta) 
      blastResults = blastSeq(refFasta, database.db)
    
    } else if (params.blastMethod == "diamond") {

      database = diamondDatabase(processPeptides.peptideFasta) 
      blastResults = diamondBlast(refFasta, database.db)
    }

    processResults = processXml(blastResults.result)

    mergeFiles = mergeeResultsFiles(processPeptides.pepResults, processResults.resultFormated)

}