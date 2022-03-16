#!/usr/bin/env nextflow
params.genera = "Brucella,Ochrobactrum,Agrobacterium"
params.download = false
params.model = "xgb"
params.k = 5
params.datadir = "$baseDir"
params.ksize = 11

nextflow.enable.dsl = 2

include { METADATA } from './workflow/metadata'
include { DOWNLOAD } from './workflow/download'
include { MAKEGRAPHS } from './workflow/makegraphs'
include { MAKEFASTA } from './workflow/makefasta'
include { PANGENOMES } from './workflow/pangenomes'
include { PRESENCE } from './workflow/presence'

workflow {
	if(params.download == true) {
		DOWNLOAD(params.datadir, params.genera)
		METADATA(params.k, DOWNLOAD.out)
	} else {
    PANGENOMES(params.datadir)
    MAKEGRAPHS(PANGENOMES.out)
    MAKEFASTA(MAKEGRAPHS.out)
    PRESENCE(MAKEFASTA.out)
    }
}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}