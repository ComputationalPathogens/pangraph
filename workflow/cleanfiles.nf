
process CLEANFILES {
	echo true
	input:
	  val(k)
	  val(datadir)

	output:
	  stdout emit: out

    script:
    """
    #!python3
    import sys
    sys.path.append("$baseDir/pyfiles/")
    import metadata
    k = "$k"
    data = "$datadir"
    out = metadata.clean_outliers(k, data)
    print(out, end = '')
    """
}