#!/usr/bin/env nextflow

//---------------------------------------
// include the RNA seq workflow
//---------------------------------------

include { epitopesBlast } from  './modules/epitopesBlast.nf'

//======================================

  if(!params.fasta) {
    throw new Exception("Missing parameter params.fasta")
  }

refFasta = Channel.fromPath(params.fasta)

  workflow { 
    epitopesBlast(refFasta)
}


