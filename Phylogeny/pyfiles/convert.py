import pandas as pd

def convert_matrix(dataloc, filenamenp = '/processed_data/features.pkl', filenamecsv = '/processed_data/cleanwcounts.csv'):
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
    datapth = dataloc + filenamenp
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
    with open('countmers.fasta', 'w') as wf:
        for s in encodedstrs:
            wf.write('>' + str(id) + '\n')
            id += 1
            wf.write(s + '\n')