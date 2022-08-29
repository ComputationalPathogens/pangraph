
process PANGENOMES {
    echo true
	input:
	  val(datadir)
	  val(dname)
          val(bifrost)
	  
	output:
	  val(datadir)

	script:
	"""
	#!python3
	import sys
	sys.path.append("$baseDir")
	from pyfiles import pangenomes
	data = "$datadir"
	out = pangenomes.pangenomes(data, "$dname", "$bifrost")
	print(out, end = '')
	"""
}
