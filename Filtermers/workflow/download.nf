

process DOWNLOAD {
	input:
	  val(datadir)
	  val(genera)
	  val(dname)

	output:
	  val(datadir)

	script:
	"""
	ncbi-genome-download --formats fasta,assembly-report  --genera "$genera" bacteria --parallel 4 -o "$datadir/samples/"
	mkdir "$datadir/samples/$dname"
	mv "$datadir/samples/refseq/bacteria"/* "$datadir/samples/$dname/"
	"""
}
