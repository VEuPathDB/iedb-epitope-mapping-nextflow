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

process makeBlastDatabase {

     publishDir "${params.blastDb}"

     container = 'veupathdb/blastsimilarity'

    input:
      path(fasta)

    output:
      path("*.fasta*")
    script:
       sample_base = fasta.getSimpleName()
       template 'makeBlastDb.bash'

}


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
/**
* mergeeResultsFiles merges the results from the exact matches and the blast into a single output

* @exactMatch output of the exact match above
* @balstOutput
*/

process mergeBlastResults {

   container = 'veupathdb/epitopemapping'

  input:
    path(blast)

  
  output:
    path("*txt")

  script:
    """
    cat ${blast} > BlastCombinedResults.txt
    """


}

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
                      
  // TODO: make the number here a param
    peptideFiles =  processPeptides.peptideFasta
                    .splitFasta(by: 1000, file:true)


  // TODO:  use database.first()?? instead of this extra param
    blastResults = blastSeq(peptideFiles, params.blastDb )
  

    processResults = processXml(blastResults.result)


  // TODO:  replace with processResults.resultFormated.collectFile(name: "mergedOutput.txt", storeDir: params.blahDir)
    blastMerge = mergeBlastResults(processResults.resultFormated.collect())
                 
    
    mergeFiles = mergeResultsFiles(processPeptides.pepResults, blastMerge)

}
