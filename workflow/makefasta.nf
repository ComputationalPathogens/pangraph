
process MAKEFASTA {
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
from pyfiles import makefasta
data = "$datadir"
out = makefasta.create(data)
print(out, end = '')
	"""
}