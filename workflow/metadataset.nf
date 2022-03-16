
process METADATASET {
    echo true
	input:
	  val(datadir)
	  val(outname)
	  val(binprefix)
	  val(numbins)

	output:
	  val(datadir)

	script:
	"""
#!python3
import sys
sys.path.append("$baseDir")
from pyfiles import metagraphs
data = "$datadir"
outname = "$outname"
binprefix = "$binprefix"
numbins = "$numbins"
out = metagraphs.build_folds(data,outname,binprefix,int(numbins))
print(out, end = '')
	"""
}