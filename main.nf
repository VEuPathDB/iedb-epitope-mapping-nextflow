#!/usr/bin/env nextflow

import nextflow.splitter.CsvSplitter

def fetchDataLink( tsv ) {
    def splitter = new CsvSplitter().options( header:true, sep:'\t' )
    def reader = new BufferedReader( new FileReader( tsv ) )
    splitter.parseHeader( reader )
    List<String> dataLink = []
    Map<String,String> row
    while( row = splitter.fetchRecord( reader ) ) {
       dataLink.add( row['dataLink'] )
    }
    return dataLink
}


//---------------------------------------
// include the RNA seq workflow
//---------------------------------------

include { iedp } from  './modules/iedbAnalysis.nf'

//======================================

if(!params.iedpSequences) {
    throw new Exception("Missing parameter params.iedpSequences")
  }
//if(!params.speciesProteome) {
 //   throw new Exception("Missing parameter params.speciesProteome")
 // }

input = fetchDataLink(params.iedpSequences)
seq = Channel.fromList(input)

  workflow {
    iedp(seq)
}