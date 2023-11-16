#!/usr/bin/env nextflow

//---------------------------------------
// include the RNA seq workflow
//---------------------------------------

include { epitopesBlast } from  './modules/epitopesBlast.nf'

//======================================

  if(!params.refFasta) {
    throw new Exception("Missing parameter params.refFasta")
  }
  if(!params.peptidesTab) {
    throw new Exception("Missing parameter params.peptidesTab")
  }
  if(!params.peptideGeneFasta) {
    throw new Exception("Missing parameter params.peptideGeneFasta")
  }

refFasta = Channel.fromPath(params.refFasta)
peptidesTab = Channel.fromPath(params.peptidesTab)
peptidesGeneFasta = Channel.fromPath(params.peptideGeneFasta)

  workflow { 
    epitopesBlast(refFasta, peptidesTab, peptidesGeneFasta)
}


