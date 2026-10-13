#!/usr/bin/env nextflow

//---------------------------------------
// include the RNA seq workflow
//---------------------------------------

include { epitopeMapping } from  './workflows/epitopeMapping.nf'

//======================================

//splitRefFasta = Channel.fromPath(params.refFasta, checkIfExists:true).splitFasta( by: params.chunkSize, file: true )

// process 10,000 epitopes at a time
//peptidesTab = Channel.fromPath(params.peptidesTab, checkIfExists: true).splitText( by: params.peptidesChunkSize, file: true )

workflow {
  if(!params.refFasta) {
    error "Missing parameter params.refFasta"
  }
  if(!params.peptidesTab) {
    error "Missing parameter params.peptidesTab"
  }
  
  if(!params.taxon) {
    error "Missing parameter params.taxon"
  }
  // if(!params.peptideMatchResults) {
  //   throw new Exception("Missing parameter params.peptideMatchResults")
  // }
  // if(!params.peptidesFilteredBySpeciesFasta) {
  //   throw new Exception("Missing parameter params.peptidesFilteredBySpeciesFasta")
  // }
  if(!params.peptideMatchResults) {
    error "Missing parameter params.peptideMatchResults"
  }
  // if(!params.chunkSize) {
  //   throw new Exception("Missing parameter params.chunkSize")
  // }
  if(!params.results) {
    error "Missing parameter params.results"
  } 

    epitopeMapping()
}


