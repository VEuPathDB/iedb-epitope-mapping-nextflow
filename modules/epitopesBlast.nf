#!/usr/bin/env nextflow
nextflow.enable.dsl=2 


process peptideSimilarity {

    publishDir "${params.results}", mode: 'copy'


    input:
    path(refFasta)
    path(pepfasta)
    path(pepTab)

    output:
    path("*Peptides.fasta"), emit: peptideFasta
    path("*PeptideGene.txt"), emit: pepResults

    script:
      template 'ProcessPeptides.bash'

}

process makeBlastDatabase {

    input:
      path(fasta)

    output:
      path("BlastDB"), emit: db
    script:
       sample_base = fasta.getSimpleName()
       template 'makeBlastDb.bash'

}


process blastSeq {

   publishDir "${params.results}/BlastOut", mode: 'copy'

   input:
    path(query)
    path(db)

   output:
   path("${sample_base}*txt")

   script:
    sample_base = query.getSimpleName()
    template 'BlastSeq.bash'

}

workflow epitopesBlast {

   take: 
    refFasta
    peptidesTab
    peptidesGeneFasta

    main:

    processPeptides = peptideSimilarity(refFasta, peptidesGeneFasta, peptidesTab)

    blastDb = makeBlastDatabase(refFasta) 
    blastResults = blastSeq(processPeptides.peptideFasta, blastDb.db)

}