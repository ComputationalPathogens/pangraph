
process KMERCOUNT {
	echo true
	input:
	  val(datadir)
	  val(dname)

	output:
	  stdout emit: dumpname

	script:
	"""
	#!python3
	import sys
	sys.path.append("$baseDir/pyfiles/")
	import kmer

	out = kmer.count_kmer("$datadir", "$dname")
	print(out, end = '')
	"""
}