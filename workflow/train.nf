
include { TRAINXGB } from './trainxgb'
include { TRAINKERAS } from './trainkeras'

workflow TRAIN {
	take:
		k
		datapth
		model

	main:
		if(model == 'xgb') {
			TRAINXGB(k, datapth)
			finalmsg = TRAINXGB.out
		} else {
			TRAINKERAS(k, datapth)
			finalmsg = TRAINKERAS.out
		}
		
    emit:
        finalmsg
		
}