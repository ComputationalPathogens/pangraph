params.datadir = "$baseDir"

include { KMERCOUNT } from './kmer'
include { BUILDMATRIX } from './buildmatrix'

workflow FEATURES {
	take:
	  datadir
	  dname

	main:
	  KMERCOUNT(datadir, dname)
	  BUILDMATRIX(KMERCOUNT.out, dname)

	emit:
	  BUILDMATRIX.out

}