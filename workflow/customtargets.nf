
process CUSTOMTARGETS {
    echo true
	input:
	  val(datadir)
	  val(targpth)

	output:
	  val(datadir)

	script:
	"""
	#!python3
import sys
sys.path.append("$baseDir")
from pyfiles import updatetargets
data = "$datadir"
targets = "$targpth"
out = updatetargets.update_targets(data,targets)
print(out, end = '')
	"""
}