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
       """
       mkdir BlastDB

       cat ${fasta} > BlastDB/${sample_base}.fa

       makeblastdb -in BlastDB/${sample_base}.fa   -title "Cookbook demo" -dbtype prot  
       """


}


process blastSeq {

   input:
    path(query)
    path(db)

   output:

   script:

   """
    blastp -query ${query} -outfmt 6 -db ${db}/*fa -out Test.txt
   """

}

workflow {

   blastDb = makeBlastDatabase(params.fasta) 
   //blastDb.db.view()
   blastResults = blastSeq(params.query, blastDb.db)

}