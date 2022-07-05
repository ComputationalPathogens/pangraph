
process PLOT {
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
	import clustermap

	out = clustermap.create("$datadir", "$dname")
	print(out, end = '')
	"""
}