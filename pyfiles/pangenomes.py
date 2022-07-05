import os
import pandas as pd
import numpy as np
from sklearn.model_selection import StratifiedKFold

def pangenomes(datadir):
    
    ####
    #Skip this step if all 5 folds have already been created, need to do this so graph splits don't get re-done every time
    ####
    if os.path.isfile(datadir + '/processed_data/fold5.gfa'):
        return datadir
    
    #Possibly make helper function for loading this in everywhere its used?
    colnames = ['id', 'assembly', 'genus', 'species', 'seqfile', 'cntfile', 'meta']
    loadpth = datadir + '/processed_data/clean.csv'
    readcsv = pd.read_csv(loadpth, names=colnames)
    paths = readcsv.seqfile.tolist()
    labels = readcsv.species.tolist()
    
    ####
    #Creating splits for both graph model and kmer model and saving to file
    ####
    kf = StratifiedKFold(n_splits=5, shuffle=True)

    test_splits = []

    for train_index, test_index in kf.split(paths, labels):
        test_splits.append(test_index)

    ###
    #Hardcoding of all 5 iterations of split combinations
    ###
    graphsplits = [np.concatenate((test_splits[0],test_splits[1],test_splits[2])),np.concatenate((test_splits[1],test_splits[2],test_splits[3])),np.concatenate((test_splits[2],test_splits[3],test_splits[4])),
                   np.concatenate((test_splits[0],test_splits[3],test_splits[4])),np.concatenate((test_splits[0],test_splits[1],test_splits[4]))]
    trainsplits = [test_splits[3],test_splits[4],test_splits[0],test_splits[1],test_splits[2]]
    testsplits = [test_splits[4],test_splits[0],test_splits[1],test_splits[2],test_splits[3]]
    
    ind = 0
    ####
    #Building pangenome graph for each fold with 60% of the dataset
    ####
    for g in graphsplits:
        ind += 1
        writepth = datadir + '/processed_data/fold' + str(ind) + 'graphsamples.txt'
        with open(writepth, 'w') as f:
            for index in g:
                f.write(datadir + paths[index] + '\n')
        cmd = '/home/liam/bifrost/bifrost/build/src/Bifrost build -r ' + writepth + ' -o ' + datadir + '/processed_data/fold' + str(ind) + ' -t 128 -c'
        os.system(cmd)
    
    
    outpth = datadir + '/processed_data/foldsplits.npy'
    np.save(outpth, [trainsplits,graphsplits,testsplits], allow_pickle=True)
    for x in range(5):
        outpth = datadir + '/processed_data/fold' + str(x+1) + 'splits.npy'
        tosave = [trainsplits[x], graphsplits[x], testsplits[x]]
        np.save(outpth, tosave)
    return datadir