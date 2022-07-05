import os
import pandas as pd

def create(datadir):
    
    ####
    #Data loader helper function?
    ####
    colnames = ['id', 'assembly', 'genus', 'species', 'seqfile', 'cntfile', 'meta']
    loadpth = datadir + '/processed_data/clean.csv'
    readcsv = pd.read_csv(loadpth, names=colnames)
    inds = readcsv.id.tolist()
    
    ####
    # Making .fasta files for the individual graph dataset, used for queries and other steps
    ####
    for ind in inds:
        openpth = datadir + '/processed_data/graphs/graph' + str(ind) + '.gfa'
        writepth = datadir + '/processed_data/fasta/graph' + str(ind) + '.fasta'
        if os.path.isfile(writepth):
            continue
        with open(openpth, 'r') as f:
            ind = 0
            with open(writepth, 'w') as w:
                for l in f:
                    temp = str.rsplit(l)
                    if temp[0] == 'S':
                        seq = '>' + str(ind)
                        w.write(seq + '\n')
                        w.write(temp[2] + '\n')
                        ind += 1
                        
    return datadir
