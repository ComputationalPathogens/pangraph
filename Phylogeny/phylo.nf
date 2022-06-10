#!/usr/bin/env nextflow
params.datadir = "$baseDir"

nextflow.enable.dsl = 2

include { INDEX } from './workflow/index'
include { COUNTS } from './workflow/counts'
include { MATRIX } from './workflow/matrix'

workflow {
    INDEX(params.datadir)
    COUNTS(INDEX.out)
    MATRIX(COUNTS.out)
}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}