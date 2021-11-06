#!/usr/bin/env nextflow


process TRAINKERAS {
	echo true
	cache false
	
	input:
	  val(k)
	  val(datapth)
	  
    output:
      stdout emit: kerasout

	script:
	"""
	#!python3
	import sys
	sys.path.append("$baseDir/pyfiles/")

	import trainmodel as tm

	data, label_encoded_y, labels_unencoded = tm.load_data("$datapth")
	final_hps, final_models, final_features, final_labels = tm.train_keras($k, data, label_encoded_y, labels_unencoded)
	tm.test_keras(final_models, final_features, final_labels)
	print("Keras Training Complete")
	"""
}
