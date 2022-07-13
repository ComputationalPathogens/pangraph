
process TRAINXGB {
    cache false
	echo true
	input:
	  val(k)
	  val(datapth)
	  val(dname)

    output:
      stdout emit: xgbout
      
	script:
	"""
	#!python3
	import sys
	sys.path.append("$baseDir/pyfiles/")
	import trainmodel
	data,label_encoded_y, labels_unencoded = trainmodel.load_data("$datapth", "$dname")
	final_models, final_features, final_labels = trainmodel.train_model($k, data, label_encoded_y, labels_unencoded, True, "$datapth", "$dname")
	print("XGB Testing Complete")
	"""
}
