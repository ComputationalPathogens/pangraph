
process BUILDMATRIX {
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
	sys.path.append("$baseDir/pyfiles/")
	import matrix
	if __name__ == '__main__':
		out = matrix.build_matrix("$datadir", "$dname")
		print(out, end = '')
	"""
}