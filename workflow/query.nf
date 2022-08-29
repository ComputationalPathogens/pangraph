
process QUERY {
    echo true
	input:
	  val(datadir)
	  val(dname)
          val(blastfrost)

	output:
	  val(datadir)

	script:
	"""
	#!python3
	import sys
	sys.path.append("$baseDir")
	from pyfiles import doqueries
	data = "$datadir"
	out = doqueries.query(data, "$dname", "$blastfrost")
	print(out, end = '')
	"""
}
