
process PLOT {
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
	import clustermap

	out = clustermap.create("$datadir")
	print(out, end = '')
	"""
}