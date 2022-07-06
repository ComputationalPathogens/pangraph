#!/usr/bin/env nextflow
params.datadir = "$baseDir"
params.virus = false
params.virusfile = ''
params.plot = false
params.dname = 'coronavirus'

nextflow.enable.dsl = 2

include { GRAPH } from './workflow/buildgraph'
include { SEPERATE } from './workflow/seperate'
include { UNITIGS } from './workflow/getunitigs'

workflow {
    if (params.virus == true) {
        SEPERATE(params.datadir, params.virusfile, params.dname)
        GRAPH(SEPERATE.out, params.dname)
    } else {
        GRAPH(params.datadir, params.dname)
    }
    UNITIGS(GRAPH.out, params.dname)
}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}
