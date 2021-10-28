import os
import pandas as pd
import numpy as np

def create(datadir):
    colnames = ['id', 'assembly', 'genus', 'species', 'seqfile', 'cntfile', 'meta']
    loadpth = datadir + '/processed_data/clean.csv'
    readcsv = pd.read_csv(loadpth, names=colnames)
    paths = readcsv.seqfile.tolist()
    inds = readcsv.id.tolist()

    for ind, s in zip(inds,paths):
        checkpth = datadir + '/processed_data/graphs/graph' + str(ind) + '.gfa'
        if os.path.isfile(checkpth):
            continue
        cmd = "/home/liam/bifrost/bifrost/build/src/Bifrost build -r "
        pth = datadir + s
        pst = ' -o ' + datadir + '/processed_data/graphs/graph' + str(ind) + ' -t 64 -c -v'
        full = cmd + pth + pst
        os.system(full)
    
    return datadir

