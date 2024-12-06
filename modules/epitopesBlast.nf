#!/usr/bin/env nextflow
nextflow.enable.dsl=2 


/**
* fetchTaxon takes an NCBI taxon ID and returns all of it child taxa.
*
* @taxonID is the NCBI taxon ID
*/

process fetchTaxon {
    container = 'veupathdb/edirect'
    
    input:
      val(taxonID)

    output:
      path("TaxaList.txt")

    script:
      template 'fetchTaxa.bash'
}

/**
*  Filter iedb tab file by a file of taxa ids.  Return the unique set of protein accessions

* @peptideTabfile: the tab file containing the peptides
* @childTaxaFile: The taxa file generated from above that contain all the child taxa for the reference being analyzed. 
*/

process peptideProteinAccessionsFilteredByTaxa {
  container = 'veupathdb/epitopemapping'

  input:
  path(peptideTabfile)
  path(childTaxaFile)

  output:
  path("filteredPeptideAccessions.txt")

  script:
  """
  peptideProteinAccessionsFilteredByTaxa.py ${peptideTabfile} ${childTaxaFile} allAccessions.tmp
  sort -u allAccessions.tmp > filteredPeptideAccessions.txt
  """

}

/*
* This step download the peptides' source proteins in which their taxa are the same as the our reference spp. 
*
* @proteinIDs: the text file containing the list of priotein IDs
*/

process fetchProtein {

    container = 'veupathdb/edirect'

    input:
      path(proteinIDs)

    output:
      path("output.fasta")


    script:
      template 'fetchProteins.bash'

}
/**
* peptideExactMatches takes a reference proteome, peptide source proteome, peptide tab file and a taxon ID to 
* identity if the reference exactly matches a peptide protein and whether the peptide match and the reference taxon is the same as peptide source taxon.
* 
* If the peptide source taxon and reference taxon are the same, the peptide is added and saved to a fasta file for the blast below. 
* @refFasta is the reference proteome fasta
* @pepfasta is the peptide proteome fasta
* @pepTab is the peptide tab file
* @taxon is the tacon of the reference
*/

process peptideExactMatches {
    container = 'veupathdb/epitopemapping'

    input:
      path(refFasta)
      path(pepProtfasta)
      path(pepTab)
      val(taxon)

    output:
      path(params.peptideMatchResults), emit: pepResults
      path(params.peptidesFilteredBySpeciesFasta), emit: peptideFasta

    script:
    """
    CheckForGene.py -r ${refFasta} \
      -l ${pepTab}  \
      -e ${pepProtfasta} \
      -t ${taxon} \
      -p $params.peptideMatchResults \
      -o $params.peptidesFilteredBySpeciesFasta
    """
}

/*
* Makes a blast database using the reference proteome.
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
* This process runs a blast using the peptide as a query and the database created during makeBlastDatabase(). Only in situation where the peptide source taxon and reference taxon is the same is the blast done
*
* @query is the peptide fasta
* @db the blast database build above 
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
* Processes the blast xml output to generate a tabular format. 
*
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
* Merges the exact match searches and the blast out to generate one file as output. 
*
* @exactMatch is the exact match file from the exact match search
* @blastOutput is the process blast output file in tabular format
* @peptideMatchBlastCombinedResults is the name of the output files
*/

process mergeResultsFiles {
    container = 'veupathdb/epitopemapping'
    
    publishDir "${params.results}", mode: 'copy'

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
    splitRefFasta
    peptidesTab

  main:

    taxonFile = fetchTaxon(params.taxon)

    peptideProteinAccessions = peptideProteinAccessionsFilteredByTaxa(peptidesTab, taxonFile)

    mergedPeptideProteins = peptideProteinAccessions.splitText( by: 500, file: true)
        | fetchProtein
        | collectFile(newLine: true)

    database = makeBlastDatabase(params.refFasta)

    // the peptideFasta output here is redundant.  it makes the same filtered epitope file for each process
    processPeptides = peptideExactMatches(splitRefFasta, mergedPeptideProteins, peptidesTab, taxonFile)

    peptideSubset = processPeptides.peptideFasta.first().splitFasta( by: params.chunkSize, file: true )

    blastResults = blastSeq(peptideSubset, database)

    processResults = processXml(blastResults.result)

    mergeBlast = processResults.resultFormated.collectFile(name: "mergedBlastOutput.txt", newLine: true)

    mergedPepResults = processPeptides.pepResults.collectFile(name: params.peptideMatchResults, newLine: true, storeDir: "${params.results}" )

    // this step merges exact match with blast results
    mergeFiles = mergeResultsFiles(mergedPepResults, mergeBlast,params.peptideMatchBlastCombinedResults)

}
