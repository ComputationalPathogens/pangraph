#!/usr/bin/env nextflow
params.datadir = "$baseDir"
params.virus = false
params.virusfile = ''
params.plot = false
params.dname = 'coronavirus'

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
        SEPERATE(params.datadir, params.virusfile, params.dname)
        INDEX(SEPERATE.out, params.dname)
    } else {
        INDEX(params.datadir, params.dname)
    }
    COUNTS(INDEX.out, params.dname)
    MATRIX(COUNTS.out, params.dname)
    if (params.plot == true) {
        PLOT(MATRIX.out, params.dname)
    }
    CONVERT(MATRIX.out, params.dname)
    IQTREE(CONVERT.out, params.dname)
    HCLUST(IQTREE.out, params.dname)
}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}
