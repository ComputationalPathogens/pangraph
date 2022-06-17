
process SEPERATE {
	input:
	  val(datadir)
      val(filename)

	output:
	  val(datadir)


	script:
	"""
	#!python3
	import sys
	import os
	sys.path.append("$baseDir/pyfiles/")
	import seperate

	out = seperate.seperate_viral("$datadir", "$filename")
	print(out, end = '')
	"""
}