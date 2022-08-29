
process MAKEFASTA {
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
from pyfiles import makefasta
data = "$datadir"
out = makefasta.create(data,"$dname")
print(out, end = '')
	"""
}
