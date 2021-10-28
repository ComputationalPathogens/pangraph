
process MAKEGRAPHS {
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
from pyfiles import creategraphs
data = "$datadir"
out = creategraphs.create(data)
print(out, end = '')
	"""
}