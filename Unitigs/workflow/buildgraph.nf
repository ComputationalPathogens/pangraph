
process GRAPH {
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
	import buildgraph

	out = buildgraph.build_index("$datadir", "$dname")
	print(out, end = '')
	"""
}