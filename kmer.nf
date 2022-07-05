#!/usr/bin/env nextflow
params.genera = "Brucella,Ochrobactrum,Agrobacterium"
params.download = true
params.model = "xgb"
params.k = 5
params.datadir = "$baseDir"
params.ksize = 31
params.features = false
params.customtargets = false
params.customtargetspath = ""

nextflow.enable.dsl = 2

include { METADATA } from './workflow/metadata'
include { FEATURES } from './workflow/features'
include { TRAIN } from './workflow/train'
include { DOWNLOAD } from './workflow/download'
include { CUSTOMTARGETS } './workflow/customtargets'

workflow {
	if(params.download == true) {
		DOWNLOAD(params.datadir, params.genera)
		METADATA(params.k, DOWNLOAD.out)
	} else {
		METADATA(params.k, params.datadir)
	}
	if(params.features == true) {
		FEATURES(METADATA.out)
		if (params.customtargets == true) {
            CUSTOMTARGETS(FEATURES.out,params.customtargetspath)
            TRAIN(params.k, CUSTOMTARGETS.out, params.model)    
        } else {
        		TRAIN(params.k, FEATURES.out, params.model)
		}
	} else {
	    if (params.customtargets == true) {
        	    CUSTOMTARGETS(METADATA.out,params.customtargetspath)
            TRAIN(params.k, CUSTOMTARGETS.out, params.model)  
	    } else {
        		TRAIN(params.k, METADATA.out, params.model)
		}
	}
}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}
