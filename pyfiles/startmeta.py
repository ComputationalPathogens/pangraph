import pandas as pd
import numpy as np
import os
import sys

def build_folds(datadir, filepth, binprefix, numbins):
    #filepth = sys.argv[1] 
    #path where bins are located where each bin has the naming convention /PATH/TO/bin.*.fa where * is the bin number
    #binprefix = sys.argv[2]
    #numbins = int(sys.argv[3])
    #datadir = '/home/liam/compare'

    ####
    #Using bifrost to create the individual graphs used in training/testing of model
    ####
    cmd1 = "/home/liam/bifrost/bifrost/build/src/Bifrost build -r "
    for x in range(numbins):
        pth = filepth + '.' + str(x+1) + '.fa'
        pst = ' -o ' + datadir + '/processed_data/' + binprefix + str(x+1) + ' -t 64 -c -v'
        full = cmd1 + pth + pst
        os.system(full)
        
        
        ####
        # Making .fasta files for the individual graph, used for queries and other steps
        ####
        openpth = datadir + '/processed_data/' + binprefix + str(x+1) +  '.gfa'
        writepth = datadir + '/processed_data/'+ binprefix + str(x+1) + '.fasta'
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
            
        for y in range(5):
            ####
            # Running BlastFrost queries to pangeome graph
            ####
            pths = '-q ' + datadir + '/processed_data/' + binprefix + str(x+1) + '.fasta ' 
            cmd = '/home/liam/BlastFrost/build/BlastFrost -g ' + datadir + '/processed_data/fold' + str(y+1) + '.gfa -f ' + datadir + '/processed_data/fold' + str(y+1) + '.bfg_colors '
            flags = '-o ' + datadir + '/processed_data' + '/fold' + str(y+1) + ' -d -t 128 -s 5.1 -v'
            full = cmd + pths + flags
            os.system(full)
    return datadir
        