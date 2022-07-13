import os
import pandas as pd

def create(datadir, dname):
    ####
    #Using bifrost to create the individual graphs used in training/testing of model
    ####
    colnames = ['id', 'assembly', 'genus', 'species', 'seqfile', 'cntfile', 'meta']
    loadpth = datadir + '/processed_data/' + str(dname) + '_clean.csv'
    readcsv = pd.read_csv(loadpth, names=colnames)
    paths = readcsv.seqfile.tolist()
    inds = readcsv.id.tolist()

    ####
    #Making regular graph dataset
    ####
    for ind, s in zip(inds,paths):
        try:
            mkpth = datadir + '/processed_data/' + str(dname) + '_graphs/'
            os.mkdir(mkpth)
        except OSError:
            pass
        checkpth = datadir + '/processed_data/' + str(dname) + '_graphs/graph' + str(ind) + '.gfa'
        if os.path.isfile(checkpth):
            continue
        cmd = "/home/liam/bifrost/bifrost/build/src/Bifrost build -r "
        pth = datadir + s
        pst = ' -o ' + datadir + '/processed_data/' + str(dname) + '_graphs/graph' + str(ind) + ' -t 64 -c -v'
        full = cmd + pth + pst
        os.system(full)
    
    return datadir

