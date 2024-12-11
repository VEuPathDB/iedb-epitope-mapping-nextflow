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
  path("otherPeptides.fasta"), emit: otherPeptidesFasta

  script:
  """
  peptideProteinAccessionsFilteredByTaxa.py ${peptideTabfile} ${childTaxaFile} allAccessions.tmp taxaFilteredPeptides.fasta smallPeptides.fasta otherPeptides.fasta
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
      path("peptidesMatchingTaxaResults.txt"), emit: peptidesMatchingTaxaAndFullSequence

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
  path(pepMatchResults)
  path(fullSeqMatchResults)
  path(pepTab)

  output:
  path("merged.gff")

  script:
  """
  cut -f 1 -d ',' $pepMatchResults $fullSeqMatchResults | sort -u >distinctPeptides
  mergeResults.py --pepMatchResults ${pepMatchResults} \
    --fullSeqMatchResults ${fullSeqMatchResults} \
    --peptideTab ${pepTab} \
    --distinctPeptides distinctPeptides \
    --nonTaxaShortPeptideCutoff ${params.nonTaxaShortPeptideCutoff} \
    --outputFile merged.gff
  """
}

process indexResults {
  container = 'biocontainers/tabix:v1.9-11-deb_cv1'
  publishDir params.results, mode: 'copy'

  input:
    path gff
    val outputFileName

  output:
    path '*.gff.gz'
    path '*.tbi'

  script:
  """
  sort -u -k1,1 -k4,4n $gff > $outputFileName
  bgzip $outputFileName
  tabix -p gff ${outputFileName}.gz
  """
}


workflow epitopeMapping {

  main:

  taxonFile = fetchTaxon(params.taxon)

  peptideProteinAccessions = peptideProteinAccessionsFilteredByTaxa(params.peptidesTab, taxonFile)

  mergedPeptideProteins = peptideProteinAccessions.accessions.splitText( by: 500, file: true)
      | fetchProtein
      | collectFile()
      | first()

  taxaAndFullSequenceResults = iedbExactMatches(params.refFasta, mergedPeptideProteins, params.peptidesTab, taxonFile)

  smallMatch = smallExactPepMatch(peptideProteinAccessions.smallPeptidesFasta, params.refFasta)
  otherMatch = exactPepMatch(peptideProteinAccessions.otherPeptidesFasta.splitFasta( by: 100000, file: true ), params.refFasta)
  taxaMatch = inexactForTaxaPeptidesPepMatch(peptideProteinAccessions.taxaFilteredPeptidesFasta.splitFasta( by: 500, file: true ), params.refFasta)

  pepMatchResults = taxaMatch.mix(otherMatch, smallMatch).collectFile()

  gff = mergeResultsFiles(pepMatchResults, taxaAndFullSequenceResults, params.peptidesTab)

  indexResults(gff, params.peptideMatchResults)
}
