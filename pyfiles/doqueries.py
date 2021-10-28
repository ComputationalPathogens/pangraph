import os
import sys
import pandas as pd
import numpy as np

def query(datadir):
    datadir = '/home/liam/compare'
    splitspth =  datadir + '/processed_data/foldsplits.npy'
    splits = np.load(splitspth, allow_pickle=True)
    for x in range(5):
        pths = ""
        for g in splits[0][x]:
            if os.path.isfile(datadir + '/processed_data/fold' + str(x+1) + '/querygraph_graph' + str(g) + '.fasta.search'):
                continue
            pths += ('-q ' + datadir + '/processed_data/fasta/graph' + str(g) + '.fasta ')
        cmd = '/home/liam/BlastFrost/build/BlastFrost -g ' + datadir + '/processed_data/fold' + str(x+1) + '.gfa -f ' + datadir + '/processed_data/fold' + str(x+1) + '.bfg_colors '
        flags = '-o ' + datadir + '/processed_data/fold' + str(x+1) + '/querygraph -d -t 128 -s 5.1 -v'
        full = cmd + pths + flags
        if pths != "":
            os.system(full)
        pths = ""
        for g in splits[2][x]:
            if os.path.isfile(datadir + '/processed_data/fold' + str(x+1) + '/querygraph_graph' + str(g) + '.fasta.search'):
                continue
            pths += ('-q ' + datadir + '/processed_data/fasta/graph' + str(g) + '.fasta ')
        full = cmd + pths + flags
        if pths != "":
            os.system(full)
    return datadir