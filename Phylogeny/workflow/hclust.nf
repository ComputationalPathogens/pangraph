
process HCLUST {
    echo true
	input:
	  val(datadir)
	  val(dname)
	  
	output:
	  val(datadir)

	script:
    """

    Rscript ${datadir}/pyfiles/hclust.R
    """
}