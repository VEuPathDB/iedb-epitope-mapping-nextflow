#!/usr/bin/env nextflow

//---------------------------------------
// include the RNA seq workflow
//---------------------------------------

include { epitopesBlast } from  './modules/epitopesBlast.nf'

//======================================

  if(!params.fasta) {
    throw new Exception("Missing parameter params.fasta")
  }

querySeq = Channel.fromPath(params.fasta)

  workflow { 
    epitopesBlast(querySeq)
}


