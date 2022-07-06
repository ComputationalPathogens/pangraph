
process UNITIGS {
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
	import unitigmatrix

	out = unitigmatrix.get_unitigs("$datadir", "$dname")
	print(out, end = '')
	"""
}