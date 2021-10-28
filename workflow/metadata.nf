
include { PROCESSFILES } from './processfiles'
include { CLEANFILES } from './cleanfiles'

workflow METADATA {
	take:
	  k
	  datadir

	main:
	  PROCESSFILES(datadir)
	  CLEANFILES(k, PROCESSFILES.out)

	emit:
	  CLEANFILES.out

}