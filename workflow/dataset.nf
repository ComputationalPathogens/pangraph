
process DATASET {
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
from pyfiles import graphdataset
data = "$datadir"
out = graphdataset.build_folds(data, "$dname")
print(out, end = '')
	"""
}