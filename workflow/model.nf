
process MODEL {
    echo true
	input:
	  val(datadir)
	  val(meta)
	  val(metapth)
	  val(dname)

	output:
	  stdout emit: filename

	script:
	"""
#!python3
import sys
sys.path.append("$baseDir")
from pyfiles import graphmodel

data = "$datadir"
meta = "$meta"
metapth = "$metapth"
out = graphmodel.build(data, meta, metapth, "$dname")
print(out, end = '')
	"""
}