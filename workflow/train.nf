
include { TRAINXGB } from './trainxgb'
include { TRAINKERAS } from './trainkeras'

workflow TRAIN {
	take:
		k
		datapth
		model
		dname

	main:
		if(model == 'xgb') {
			TRAINXGB(k, datapth, dname)
			finalmsg = TRAINXGB.out
		} else {
			TRAINKERAS(k, datapth)
			finalmsg = TRAINKERAS.out
		}
		
    emit:
        finalmsg
		
}