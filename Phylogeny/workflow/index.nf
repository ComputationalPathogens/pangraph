
process INDEX {
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
	import indexfiles

	out = indexfiles.build_index("$datadir")
	print(out, end = '')
	"""
}