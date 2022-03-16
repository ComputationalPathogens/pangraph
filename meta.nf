#!/usr/bin/env nextflow
params.binpath = ""
params.binprefix = "metabin"
params.datadir = "$baseDir"
params.outname = "metabin"
params.numbins = 1

nextflow.enable.dsl = 2


include { METAPROCESS } from './workflow/metaprocess'
include { METADATASET } from './workflow/metadataset'

workflow {
	METAPROCESS(params.datadir, params.binpath, params.binprefix, params.numbins)
    METADATASET(METAPROCESS.out, params.outname, params.binprefix, params.numbins)
}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}