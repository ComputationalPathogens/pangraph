
process PANGENOMES {
    echo true
	input:
	  val(datadir)

	output:
	  val(datadir)

	script:
	"""
	#!python3
	import sys
	sys.path.append("$baseDir")
	from pyfiles import pangenomes
	data = "$datadir"
	out = pangenomes.pangenomes(data)
	print(out, end = '')
	"""
}