
process COUNTS {
	input:
	  val(datadir)

	output:
	  val(datadir)


	script:
	"""
	#!python3
	import sys
	import os
	sys.path.append("$baseDir/pyfiles/")
	import counts

	out = counts.count_kmer("$datadir")
	print(out, end = '')
	"""
}