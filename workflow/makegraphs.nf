
process MAKEGRAPHS {
    echo true
	input:
	  val(datadir)
	  val(dname)

	output:
	  val(datadir)

	script:
	"""
#!python3
import sys
sys.path.append("$baseDir")
from pyfiles import creategraphs
data = "$datadir"
out = creategraphs.create(data, "$dname")
print(out, end = '')
	"""
}