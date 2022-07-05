
process CONVERT {
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
	import convert

	out = convert.convert_matrix("$datadir", "$dname")
	print(out, end = '')
	"""
}