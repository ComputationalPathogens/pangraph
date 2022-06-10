
process MATRIX {
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
	import matrix

	out = matrix.build_matrix("$datadir")
	print(out, end = '')
	"""
}