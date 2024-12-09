#!/usr/bin/env nextflow

//---------------------------------------
// include the RNA seq workflow
//---------------------------------------

include { epitopeMapping } from  './workflows/epitopeMapping.nf'

//======================================

  if(!params.refFasta) {
    throw new Exception("Missing parameter params.refFasta")
  }
  if(!params.peptidesTab) {
    throw new Exception("Missing parameter params.peptidesTab")
  }
  
  if(!params.taxon) {
    throw new Exception("Missing parameter params.taxon")
  }
  // if(!params.peptideMatchResults) {
  //   throw new Exception("Missing parameter params.peptideMatchResults")
  // }
  // if(!params.peptidesFilteredBySpeciesFasta) {
  //   throw new Exception("Missing parameter params.peptidesFilteredBySpeciesFasta")
  // }
  if(!params.peptideMatchBlastCombinedResults) {
    throw new Exception("Missing parameter params.peptideMatchBlastCombinedResults")
  }
  if(!params.chunkSize) {
    throw new Exception("Missing parameter params.chunkSize")
  } 
  if(!params.results) {
    throw new Exception("Missing parameter params.results")
  } 

//splitRefFasta = Channel.fromPath(params.refFasta, checkIfExists:true).splitFasta( by: params.chunkSize, file: true )

// process 10,000 epitopes at a time
//peptidesTab = Channel.fromPath(params.peptidesTab, checkIfExists: true).splitText( by: params.peptidesChunkSize, file: true )

workflow {
    epitopeMapping()
}


