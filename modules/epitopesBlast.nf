#!/usr/bin/env nextflow
nextflow.enable.dsl=2 


process peptideSimilarity {
    input:
    path(refFasta)
    path(pepfasta)
    path(pepTab)

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

   //blastDb = makeBlastDatabase(fasta) 
   //blastResults = blastSeq(params.query, blastDb.db)

}