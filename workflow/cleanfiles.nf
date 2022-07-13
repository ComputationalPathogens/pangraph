
process CLEANFILES {
	echo true
	input:
	  val(k)
	  val(datadir)
	  val(dname)

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
    out = metadata.clean_outliers(k, data, "$dname")
    print(out, end = '')
    """
}