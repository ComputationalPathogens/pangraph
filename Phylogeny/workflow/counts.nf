
process COUNTS {
	input:
	  val(datadir)
	  val(dname)

	output:
	  val(datadir)


	script:
	"""
	#!python3
	import sys
	import os
	sys.path.append("$baseDir/pyfiles/")
	import counts

	out = counts.count_kmer("$datadir", 31, "$dname")
	print(out, end = '')
	"""
}