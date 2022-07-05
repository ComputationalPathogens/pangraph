#!/usr/bin/env nextflow
params.genera = "Bacillus anthracis"
params.download = false
params.model = "xgb"
params.k = 5
params.datadir = "$baseDir"
params.ksize = 11
params.meta = false
params.metapth = ""
params.customtargets = false
params.customtargetspath = ""

nextflow.enable.dsl = 2

include { METADATA } from './workflow/metadata'
include { DOWNLOAD } from './workflow/download'
include { MAKEGRAPHS } from './workflow/makegraphs'
include { MAKEFASTA } from './workflow/makefasta'
include { QUERY } from './workflow/query'
include { PANGENOMES } from './workflow/pangenomes'
include { DATASET } from './workflow/dataset'
include { MODEL } from './workflow/model'
include { CUSTOMTARGETS } './workflow/customtargets'

workflow {
	if(params.download == true) {
		DOWNLOAD(params.datadir, params.genera)
		METADATA(params.k, DOWNLOAD.out)
        
	} else {
    METADATA(params.k, params.datadir)
    if (params.customtargets == true) {
        CUSTOMTARGETS(METADATA.out,params.customtargetspath)
    }
    PANGENOMES(METADATA.out)
    MAKEGRAPHS(PANGENOMES.out)
    MAKEFASTA(MAKEGRAPHS.out)
    QUERY(MAKEFASTA.out)
    DATASET(QUERY.out)
    MODEL(DATASET.out, params.meta, params.metapth)
    }
}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}
