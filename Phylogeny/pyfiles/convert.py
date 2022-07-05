import pandas as pd
import os

def convert_matrix(dataloc, dname):
    """
    Parameters
    ----------
    dataloc : Base directory of nextflow execution

    Returns
    -------
    data : Feature matrix of [NumSamples]*[NumFeatures] shape
    labels_encoded : Labels corresponding to feature matrix [Numsamples] length
    label_encoder.classes_ : the unencoded classes the model is being trained on
    
    """
    if os.path.isfile(dataloc + '/processed_data/' + str(dname) + '_kmercounts.fasta'):
        return dataloc
    datapth = dataloc + '/processed_data/' + str(dname) + '_featuresfiltered.pkl'
    data = pd.read_pickle(datapth)
    matrix = data.values.tolist()
    encodedstrs = []
    for ls in matrix:
        tempstr = ''
        for n in ls:
            if n > 15:
                tempstr += 'f'
            else:
                tempstr += hex(n)[2]
        encodedstrs.append(tempstr)
    id = 0
    with open(dataloc + '/processed_data/' + str(dname) + '_kmercounts.fasta', 'w') as wf:
        for s in encodedstrs:
            wf.write('>' + str(id) + '\n')
            id += 1
            wf.write(s + '\n')
    return dataloc
