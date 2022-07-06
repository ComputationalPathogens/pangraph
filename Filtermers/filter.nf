#!/usr/bin/env nextflow
params.datadir = "$baseDir"
params.virus = false
params.virusfile = ''
params.dname = 'coronavirus'
params.download = false
params.genera = ""

nextflow.enable.dsl = 2

include { INDEX } from './workflow/index'
include { COUNTS } from './workflow/counts'
include { MATRIX } from './workflow/matrix'
include { SEPERATE } from './workflow/seperate'
include { DOWNLOAD } from './workflow/download'

workflow {
    if (params.virus == true) {
        SEPERATE(params.datadir, params.virusfile, params.dname)
        INDEX(SEPERATE.out, params.dname)
    } else if (params.download == true) {
        DOWNLOAD(params.datadir, params.genera, params.dname)
        INDEX(DOWNLOAD.out, params.dname)
    } else {
        INDEX(params.datadir, params.dname)
    }
    COUNTS(INDEX.out, params.dname)
    MATRIX(COUNTS.out, params.dname)
}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}
