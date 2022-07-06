
process SEPERATE {
	input:
	  val(datadir)
      val(filename)
      val(dname)

	output:
	  val(datadir)


	script:
	"""
	#!python3
	import sys
	import os
	sys.path.append("$baseDir/pyfiles/")
	import seperate

	out = seperate.seperate_viral("$datadir", "$filename", "$dname")
	print(out, end = '')
	"""
}