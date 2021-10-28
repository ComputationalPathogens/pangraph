

process DOWNLOAD {
	input:
	  val(datadir)
	  val(genera)

	output:
	  val(datadir)

	script:
	"""
	ncbi-genome-download --formats fasta,assembly-report  --genera "$genera" bacteria --parallel 4 -o "$datadir"
	"""
}