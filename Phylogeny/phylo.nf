#!/usr/bin/env nextflow
params.datadir = "$baseDir"
params.virus = false
params.virusfile = ''
params.plot = false

nextflow.enable.dsl = 2

include { INDEX } from './workflow/index'
include { COUNTS } from './workflow/counts'
include { MATRIX } from './workflow/matrix'
include { SEPERATE } from './workflow/seperate'
include { CONVERT } from './workflow/convert'
include { PLOT } from './workflow/plot'
include { IQTREE } from './workflow/iqtree'
include { HCLUST } from './workflow/hclust'

workflow {
    if (params.virus == true) {
        SEPERATE(params.datadir, params.virusfile)
        INDEX(SEPERATE.out)
    } else {
        INDEX(params.datadir)
    }
    COUNTS(INDEX.out)
    MATRIX(COUNTS.out)
    if (params.plot == true) {
        PLOT(MATRIX.out)
    }
    CONVERT(MATRIX.out)
    IQTREE(CONVERT.out)
    HCLUST(IQTREE.out)
}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}
