import os
import pandas as pd

def create(datadir, dname):
    
    ####
    #Data loader helper function?
    ####
    colnames = ['id', 'assembly', 'genus', 'species', 'seqfile', 'cntfile', 'meta']
    loadpth = datadir + '/processed_data/' + str(dname) + '_clean.csv'
    readcsv = pd.read_csv(loadpth, names=colnames)
    inds = readcsv.id.tolist()
    
    ####
    # Making .fasta files for the individual graph dataset, used for queries and other steps
    ####
    try:
        mkpth = datadir + '/processed_data/' + str(dname) + '_fasta/'
        os.mkdir(mkpth)
    except OSError:
        pass
    for ind in inds:
        openpth = datadir + '/processed_data/' + str(dname) + '_graphs/graph' + str(ind) + '.gfa'
        writepth = datadir + '/processed_data/' + str(dname) + '_fasta/graph' + str(ind) + '.fasta'
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
