
process HCLUST {
    echo true
	input:
	  val(datadir)

	output:
	  val(datadir)

	script:
    """

    Rscript ${datadir}/pyfiles/hclust.R
    """
}