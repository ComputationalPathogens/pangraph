
process PRESENCE {
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
from pyfiles import unitigpresence
data = "$datadir"
out = unitigpresence.get_presence(data)
print(out, end = '')
	"""
}
