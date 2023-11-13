#!/usr/bin/env nextflow
nextflow.enable.dsl=2 

params.fasta = "$projectDir/data/query/query.fa"
params.query = "$projectDir/Results/iedbFasta/iedbEpitpes.fa"

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

   input:
    path(query)
    path(db)

   output:

   script:
    template 'BlastSeq.sh'

   //"""
   // blastp -query ${query} -outfmt 6 -db ${db}/*fa -out Test.txt
   //"""

}

workflow epitopesBlast {

   take: 
    querySeq

    main:

   blastDb = makeBlastDatabase(params.fasta) 
   blastResults = blastSeq(params.query, blastDb.db)

}