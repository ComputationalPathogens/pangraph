
process INDEX {
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
	import indexfiles

	out = indexfiles.build_index("$datadir", "$dname")
	print(out, end = '')
	"""
}