
process MATRIX {
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
	import matrix

	out = matrix.build_matrix("$datadir", "$dname")
	print(out, end = '')
	"""
}