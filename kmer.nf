#!/usr/bin/env nextflow
params.genera = "Brucella,Ochrobactrum,Agrobacterium"
params.download = false
params.model = "xgb"
params.k = 5
params.datadir = "$baseDir"
params.ksize = 31
params.features = false
params.customtargets = false
params.customtargetspath = ""
params.dname = ""

nextflow.enable.dsl = 2

include { METADATA } from './workflow/metadata'
include { FEATURES } from './workflow/features'
include { TRAIN } from './workflow/train'
include { DOWNLOAD } from './workflow/download'
include { CUSTOMTARGETS } from './workflow/customtargets'

workflow {
	if(params.download == true) {
		DOWNLOAD(params.datadir, params.genera, params.dname)
		METADATA(params.k, DOWNLOAD.out, params.dname)
	} else {
		METADATA(params.k, params.datadir, params.dname)
	}
	if(params.features == true) {
		FEATURES(METADATA.out, params.dname)
		if (params.customtargets == true) {
            CUSTOMTARGETS(FEATURES.out,params.customtargetspath, params.dname)
            TRAIN(params.k, CUSTOMTARGETS.out, params.model, params.dname)    
        } else {
        		TRAIN(params.k, FEATURES.out, params.model, params.dname)
		}
	} else {
	    if (params.customtargets == true) {
        	    CUSTOMTARGETS(METADATA.out,params.customtargetspath, params.dname)
            TRAIN(params.k, CUSTOMTARGETS.out, params.model, params.dname)  
	    } else {
        		TRAIN(params.k, METADATA.out, params.model, params.dname)
		}
	}
}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}
