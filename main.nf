#!/usr/bin/env nextflow

//---------------------------------------
// include the RNA seq workflow
//---------------------------------------

include { iedp } from  './modules/iedbAnalysis.nf'

//======================================

  if(!params.iedpUrl) {
    throw new Exception("Missing parameter params.iedpUrl")
  }

url = Channel.fromPath(params.iedpUrl)

  workflow {
    iedp(url)
}