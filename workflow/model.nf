
process MODEL {
    echo true
	input:
	  val(datadir)

	output:
	  stdout emit: filename

	script:
	"""
#!python3
import sys
sys.path.append("$baseDir")
from pyfiles import graphmodel

data = "$datadir"
out = graphmodel.build(data)
print(out, end = '')
	"""
}