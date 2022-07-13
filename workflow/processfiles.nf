
process PROCESSFILES {
	echo true
	input:
	  val(datadir)
	  val(dname)

	output:
	  stdout emit: out


	script:
	"""
	#!python3
	import sys
	import os
	sys.path.append("$baseDir/pyfiles/")
	import metadata

	out = metadata.build_metadata("$datadir", "$dname")
	print(out, end = '')
	"""
}