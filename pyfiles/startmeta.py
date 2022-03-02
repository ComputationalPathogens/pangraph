import pandas as pd
import numpy as np
import os
import sys

def build_folds():
    filepth = sys.argv[1]
    outname = sys.argv[2]
    datadir = '/home/liam/compare'

    ####
    #Using bifrost to create the individual graphs used in training/testing of model
    ####
    colnames = ['id', 'assembly', 'genus', 'species', 'seqfile', 'cntfile', 'meta']
    loadpth = datadir + '/processed_data/clean.csv'
    readcsv = pd.read_csv(loadpth, names=colnames)
    paths = readcsv.seqfile.tolist()
    inds = readcsv.id.tolist()

    ####
    #Making regular graph dataset
    ####
    cmd = "/home/liam/bifrost/bifrost/build/src/Bifrost build -r "
    pth = filepth
    pst = ' -o ' + datadir + '/processed_data/' + outname + ' -t 64 -c -v'
    full = cmd + pth + pst
    os.system(full)
    
    
    ####
    # Making .fasta files for the individual graph dataset, used for queries and other steps
    ####
    openpth = datadir + '/processed_data/' + outname +  '.gfa'
    writepth = datadir + '/processed_data/'+ outname + '.fasta'
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
        
    for x in range(5):
        ####
        # Regular graph dataset
        ####
        pths = '-q ' + datadir + '/processed_data/' + outname + '.fasta ' 
        cmd = '/home/liam/BlastFrost/build/BlastFrost -g ' + datadir + '/processed_data/fold' + str(x+1) + '.gfa -f ' + datadir + '/processed_data/fold' + str(x+1) + '.bfg_colors '
        flags = '-o ' + datadir + '/processed_data' + '/fold' + str(x+1) + ' -d -t 128 -s 5.1 -v'
        full = cmd + pths + flags
        os.system(full)
        
build_folds()