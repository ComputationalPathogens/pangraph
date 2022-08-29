#!/usr/bin/env nextflow
params.genera = "Brucella"
params.download = false
params.k = 5
params.datadir = "$baseDir"
params.meta = false
params.metapth = ""
params.customtargets = false
params.customtargetspath = ""
params.dname = ""
params.bifrost = ""
params.blastfrost = ""

nextflow.enable.dsl = 2

include { METADATA } from './workflow/metadata'
include { DOWNLOAD } from './workflow/download'
include { MAKEGRAPHS } from './workflow/makegraphs'
include { MAKEFASTA } from './workflow/makefasta'
include { QUERY } from './workflow/query'
include { PANGENOMES } from './workflow/pangenomes'
include { DATASET } from './workflow/dataset'
include { MODEL } from './workflow/model'
include { CUSTOMTARGETS } from './workflow/customtargets'

workflow {
	if(params.download == true) {
		DOWNLOAD(params.datadir, params.genera, params.dname)
		METADATA(params.k, DOWNLOAD.out, params.dname)
        
	} else {
    METADATA(params.k, params.datadir, params.dname)
    if (params.customtargets == true) {
        CUSTOMTARGETS(METADATA.out,params.customtargetspath, params.dname)
    }
    PANGENOMES(METADATA.out, params.dname, params.bifrost)
    MAKEGRAPHS(PANGENOMES.out, params.dname)
    MAKEFASTA(MAKEGRAPHS.out, params.dname)
    QUERY(MAKEFASTA.out, params.dname,params.blastfrost)
    DATASET(QUERY.out, params.dname)
    MODEL(DATASET.out, params.meta, params.metapth, params.dname)
    }
}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}
