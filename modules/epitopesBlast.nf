#!/usr/bin/env nextflow
nextflow.enable.dsl=2 


/**
* fetchTaxon takes in an NCBI taxon ID and retun all of it child taxa.
*
* @taxonID is the NCBI taxon ID
*/

process fetchTaxon {
    
    container = 'veupathdb/edirect'
    input:
      val(taxonID)

    output:
      path("TaxaList.txt"), emit: taxaFile

    script:
      template 'fetchTaxa.bash'


}

/**
* The process below take in reference proteome, peptide source proteome, peptide tab file and a taxon ID to 
* identity if a the reference exactly matches a peptide protein and weather the peptide match and the referece taxon is the same as peptide source taxon.
* 
* If the peptide source taxon and reference taxon are the same, the peptide is added and save into a fasta file for the blast below. 
* @refFasta is the reference proteome fasta
* @pepfasta is the peptide proteome fasta
* @pepTab is the peptide tab file
* @taxon is the tacon of the reference
*/

process peptideExactMatches {

    container = 'veupathdb/epitopemapping'

    publishDir "${params.results}", mode: 'copy', pattern: "*txt"


    input:
      path(refFasta)
      path(pepProtfasta)
      path(pepTab)
      val(taxon)
      val(peptideMatchResultsOutput)
      val(peptidesFilteredBySpeciesFastaOutput)

    output:
      path(peptidesFilteredBySpeciesFastaOutput), emit: peptideFasta
      path(peptideMatchResultsOutput), emit: pepResults


    script:
      template 'ProcessPeptides.bash'

}

/*
* Make a blast database using the reference proteome.
*
* @fasta is the input proteome.
*/
process makeBlastDatabase {

     container = 'veupathdb/blastsimilarity'

    input:
      path(fasta)

    output:
      path("db")
    script:
       sample_base = fasta.getSimpleName()
       template 'makeBlastDb.bash'

}

/*
* This process does the blast using the peptide as query and the database created above.
*
* @query; is the peptide fasta
* @db; the blast database build above 
*/

process blastSeq {

    container = 'veupathdb/blastsimilarity'

    input:
      path(query)
      path(db)

    output:
      path("${sample_base}*xml"), emit: result

    script:
      sample_base = query.getName()
      template 'BlastSeq.bash'

}

/*
* Process the blast xml output to generate a tabular format/

* @xml is the blast xml output file to be used in subsequent steps below

*/

process processXml {

     container = 'veupathdb/epitopemapping'

    input:
      path(xml)

    output:
      path("*txt"), emit: resultFormated


    script:
      sample_base = xml.getName()
      template 'processBlastXml.bash'
}

/*

* This process the merges the exact match searches and the blast out to generate one file as output. 

* @exactMatch is the exact match file from the exact match search
* @balstOutput is the process blast output file in tabular format
* @peptideMatchBlastCombiedResults is the name of the output files

*/
process mergeResultsFiles {

    container = 'veupathdb/epitopemapping'
    
    publishDir "${params.results}", mode: 'copy', overwrite: false


    input:
      path(exactMatch)
      path(balstOutput)
      val(peptideMatchBlastCombiedResults)
    
  

    output:
      path(params.peptideMatchBlastCombinedResults)

    script:
      template 'mergeFiles.bash'

}

workflow epitopesBlast {

   take: 
    refFasta
    peptidesTab
    peptidesGeneFasta

    main:

    taxonFile = fetchTaxon(params.taxon)

    database = makeBlastDatabase(refFasta)


    processPeptides = peptideExactMatches(refFasta, peptidesGeneFasta, peptidesTab, taxonFile.taxaFile, params.peptideMatchResults, params.peptidesFilteredBySpeciesFasta)
                      

    peptideFiles =  processPeptides.peptideFasta
                    .splitFasta(by: params.chuckSize, file:true)


    blastResults = blastSeq(peptideFiles, database.first())
  

    processResults = processXml(blastResults.result)

    mergeBlast = processResults.resultFormated.collectFile(name: "mergedOutput.txt", newLine: true )
    
    mergeFiles = mergeResultsFiles(processPeptides.pepResults, mergeBlast,params.peptideMatchBlastCombinedResults)

}
