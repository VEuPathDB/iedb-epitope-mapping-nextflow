#!/usr/bin/env nextflow
nextflow.enable.dsl=2


include { pepMatch as smallExactPepMatch } from  '../subworkflows/local/pepmatch.nf'
include { pepMatch as exactPepMatch } from  '../subworkflows/local/pepmatch.nf'
include { pepMatch as inexactForTaxaPeptidesPepMatch } from  '../subworkflows/local/pepmatch.nf'

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
    """
    esearch -db taxonomy -query "txid${taxonID}[Subtree]" \
        | efetch -format xml \
        | xtract -pattern TaxaSet -block "*/Taxon" \
                 -tab "\n" -element TaxId > TaxaList.txt
    """
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
  path("filteredPeptideAccessions.txt"), emit: accessions
  path("taxaFilteredPeptides.fasta"), emit: taxaFilteredPeptidesFasta

  path("smallPeptides.fasta"), emit: smallPeptidesFasta
  path("notSmallPeptides.fasta"), emit: notSmallPeptidesFasta

  script:
  """
  peptideProteinAccessionsFilteredByTaxa.py ${peptideTabfile} ${childTaxaFile} allAccessions.tmp taxaFilteredPeptides.fasta smallPeptides.fasta notSmallPeptides.fasta
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
* @refFasta is the reference proteome fasta
* @pepfasta is the peptide proteome fasta
* @pepTab is the peptide tab file
* @taxon is the tacon of the reference
*/

process iedbExactMatches {
    container = 'veupathdb/epitopemapping'

    input:
      path(refFasta)
      path(pepProtfasta)
      path(pepTab)
      path(taxaFile)

    output:
      path("peptidesMatchingTaxaResults.txt"), emit: peptidesMatchingTaxaResults

    script:
    """
    processIedbExactMatches.py --refProteome ${refFasta} \
      --epitopetab ${pepTab}  \
      --epitopeProtein ${pepProtfasta} \
      --refTaxa ${taxaFile} \
      --peptideMatchOutput peptidesMatchingTaxaResults.txt \
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
      path("${sample_base}*xml")

    script:
      sample_base = query.getName()
      template 'BlastSeq.bash'
}

/*
* Processes the blast xml output to generate a tabular format. 
*
* @xml is the blast xml output file to be used in subsequent steps below
*/

process parseBlastXml {
    container = 'veupathdb/epitopemapping'

    input:
      path(xml)

    output:
      path("parsedBlast.txt")

    script:
    """
    parseXML.pl --inputXml ${xml} \
                --outputFile parsedBlast.txt
    """
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
    """"
    mergeBlastAndExactMatch.pl --exactMatchFile ${exactMatch} \
                               --blastFile ${balstOutput} \
                               --outputFile ${peptideMatchBlastCombiedResults}
    """
}

workflow epitopeMapping {

  // take:
  //   splitPeptidesTab

  main:
//  database = makeBlastDatabase(params.refFasta)
  taxonFile = fetchTaxon(params.taxon)

  peptideProteinAccessions = peptideProteinAccessionsFilteredByTaxa(params.peptidesTab, taxonFile)

  // mergedPeptideProteins = peptideProteinAccessions.accessions.splitText( by: 500, file: true)
  //     | fetchProtein
  //     | collectFile()
  //     | first()

//  processPeptides = iedbExactMatches(params.refFasta, mergedPeptideProteins, params.peptidesTab, taxonFile)

  smallExactPepMatch(peptideProteinAccessions.smallPeptidesFasta, params.refFasta)
  exactPepMatch(peptideProteinAccessions.notSmallPeptidesFasta.splitFasta( by: 100000, file: true ), params.refFasta)
  inexactForTaxaPeptidesPepMatch(peptideProteinAccessions.taxaFilteredPeptidesFasta.splitFasta( by: 500, file: true ), params.refFasta)

  // mergedPepResults = processPeptides.pepResults.collectFile()

    // peptideSubset = peptideProteinAccessions.peptides.splitFasta( by: params.chunkSize, file: true )
    // mergedBlastResults = blastSeq(peptideSubset, database)
    //      | parseBlastXml
    //      | collectFile(newLine: true)


    // // this step merges exact match with blast results
    // mergeFiles = mergeResultsFiles(mergedPepResults, mergedBlastResults, params.peptideMatchBlastCombinedResults)

}
