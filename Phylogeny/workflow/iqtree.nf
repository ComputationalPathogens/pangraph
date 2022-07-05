
process IQTREE {
	input:
	  val(datadir)

	output:
	  val(datadir)


	script:
	"""
    /home/liam/iqtree-1.6.12-Linux/bin/iqtree -m GTR2+FO+ASC+R9 -s ${datadir}/processed_data/kmercounts.fasta -nt AUTO -mem 150G
	"""
}
