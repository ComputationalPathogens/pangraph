
include { PROCESSFILES } from './processfiles'
include { CLEANFILES } from './cleanfiles'

workflow METADATA {
	take:
	  k
	  datadir
	  dname

	main:
	  PROCESSFILES(datadir, dname)
	  CLEANFILES(k, PROCESSFILES.out, dname)

	emit:
	  CLEANFILES.out

}