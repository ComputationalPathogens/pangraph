
process METAPROCESS {
    echo true
	input:
	  val(datadir)
	  val(binpath)
	  val(binprefix)
	  val(numbins)

	output:
	  val(datadir)

	script:
	"""
#!python3
import sys
sys.path.append("$baseDir")
from pyfiles import startmeta
data = "$datadir"
binpath = "$binpath"
binprefix = "$binprefix"
numbins = "$numbins"
out = startmeta.build_folds(data,binpath,binprefix,int(numbins))
print(out, end = '')
	"""
}