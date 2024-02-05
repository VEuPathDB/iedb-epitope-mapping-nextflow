#!/usr/bin/env nextflow
nextflow.enable.dsl=2 


/**
* fetchTaxon takes in an NCBI taxon ID and retun all of it child taxa.
*
* @taxonID is the NCBI taxon ID
*/

process fetchTaxon {
    //FIXME:  Update the containter here to the veupath one
    container = 'ncbi/edirect:latest'
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

    //container = 'veupathdb/epitopemapping'
    container = 'epitopemapping'


    publishDir "${params.results}", mode: 'copy', pattern: "*txt"


    input:
      path(refFasta)
      path(pepfasta)
      path(pepTab)
      val(taxon)
      val(peptideMatchResults)
      val(peptidesFilteredBySpeciesFasta)

    output:
      path("*fasta"), emit: peptideFasta
      path("*txt"), emit: pepResults


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
    

    //publishDir "${params.results}/BlastOut", mode: 'copy'

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

     container = 'epitopemapping'

    //publishDir "${params.results}/BlastOut", mode: 'copy'

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

    container = 'epitopemapping'
    
    publishDir "${params.results}", mode: 'copy'

    input:
      path(exactMatch)
      path(balstOutput)
      val(peptideMatchBlastCombiedResults)
    
  

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

   taxonFile = fetchTaxon(params.taxon)

   database = makeBlastDatabase(refFasta) 

    processPeptides = peptideExactMatches(refFasta, peptidesGeneFasta, peptidesTab, taxonFile.taxaFile, params.peptideMatchResults, params.peptidesFilteredBySpeciesFasta)
                      
                      
    peptideFiless =  processPeptides.peptideFasta
                    .splitFasta(by: 1000, file:true)

    blastResults = blastSeq(peptideFiless, params.blastDb )
  

    processResults = processXml(blastResults.result)


    blastMerge = mergeBlastResults(processResults.resultFormated.collect())
                 
    
    mergeFiles = mergeResultsFiles(processPeptides.pepResults, blastMerge, params.peptideMatchBlastCombinedResults)

}
