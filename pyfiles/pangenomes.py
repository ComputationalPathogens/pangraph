import os
import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold

def pangenomes(datadir):
    checkpth = datadir + '/processed_data/fold5.gfa'
    if os.path.isfile(checkpth):
        return datadir
    colnames = ['id', 'assembly', 'genus', 'species', 'seqfile', 'cntfile', 'meta']
    loadpth = datadir + '/processed_data/clean.csv'
    readcsv = pd.read_csv(loadpth, names=colnames)
    paths = readcsv.seqfile.tolist()
    labels = readcsv.species.tolist()
    kf = StratifiedKFold(n_splits=5, shuffle=True)
    train_splits = []
    graph_splits = []
    test_splits = []
    for train_index, test_index in kf.split(paths, labels):
        test_splits.append(test_index)
        gsplit = len(train_index) // 4
        train = train_index[:gsplit]
        graph = train_index[gsplit:]
        train_splits.append(train)
        graph_splits.append(graph)
    ind = 0
    for g in graph_splits:
        ind += 1
        writepth = datadir + '/processed_data/fold' + str(ind) + 'graphsamples.txt'
        with open(writepth, 'w') as f:
            for index in g:
                f.write(datadir + paths[index] + '\n')
        cmd = '/home/liam/bifrost/bifrost/build/src/Bifrost build -r ' + writepth + ' -o ' + datadir + '/processed_data/fold' + str(ind) + ' -t 128 -c'
    
        os.system(cmd)
    outpth = datadir + '/processed_data/foldsplits.npy'
    np.save(outpth, [train_splits,graph_splits,test_splits], allow_pickle=True)
    for x in range(5):
        outpth = datadir + '/processed_data/fold' + str(x+1) + 'splits.npy'
        tosave = [train_splits[x], graph_splits[x], test_splits[x]]
    
        np.save(outpth, tosave)
        
            
        
        return datadir