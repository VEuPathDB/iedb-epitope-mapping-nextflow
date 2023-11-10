#!/usr/bin/env nextflow

//---------------------------------------
// include the RNA seq workflow
//---------------------------------------

include { iedbEpitopesDownload } from  './modules/iedbEpitopesDownload.nf'

//======================================

  if(!params.iedbUrl) {
    throw new Exception("Missing parameter params.iedpUrl")
  }

url = Channel.fromPath(params.iedbUrl)

  workflow {
    iedbEpitopesDownload(url)
}