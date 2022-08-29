import os
import sys
import pandas as pd
import numpy as np
import subprocess

def query(datadir, dname):
    ####
    # Building the commands for querying the pangenome graphs with the .fasta of indivudal graphs
    # Queries are done using BlastFrost (longest part of strategy, look for optimizations here)
    # Queries are done for each fold with the corresponding train/test set being queried against that folds pangenome graph
    ####
    splitspth =  datadir + '/processed_data/' + str(dname) + '_foldsplits.npy'
    splits = np.load(splitspth, allow_pickle=True)
    for x in range(5):
        ####
        # Regular graph dataset
        ####
        try:
            mkpth = datadir + '/processed_data/' + str(dname) + '_fold' + str(x+1) + '/'
            os.mkdir(mkpth)
        except:
            pass
        cmd = '/home/liam/BlastFrost/build/BlastFrost -g ' + datadir + '/processed_data/' + str(dname) + '_fold' + str(x+1) + '.gfa -f ' + datadir + '/processed_data/' + str(dname) + '_fold' + str(x+1) + '.bfg_colors '
        flags = '-o ' + datadir + '/processed_data/' + str(dname) + '_fold' + str(x+1) + '/querygraph -d -t 128 -s 5.1 -v'
        pths = ""
        filecnt = 0
        for g in splits[0][x]:
            if os.path.isfile(datadir + '/processed_data/' + str(dname) + '_fold' + str(x+1) + '/querygraph_graph' + str(g) + '.fasta.search'):
                continue
            pths += ('-q ' + datadir + '/processed_data/' + str(dname) + '_fasta/graph' + str(g) + '.fasta ')
            filecnt += 1
            if filecnt == 500:
                full = cmd + pths + flags
                if pths != "":
                    subprocess.call(full,shell=True)
                pths = ""
                filecnt = 0
        full = cmd + pths + flags
        if pths != "":
            subprocess.call(full,shell=True)
        pths = ""
        filecnt = 0
        for g in splits[2][x]:
            if os.path.isfile(datadir + '/processed_data/' + str(dname) + '_fold' + str(x+1) + '/querygraph_graph' + str(g) + '.fasta.search'):
                continue
            pths += ('-q ' + datadir + '/processed_data/' + str(dname) + '_fasta/graph' + str(g) + '.fasta ')
            filecnt += 1
            if filecnt == 500:
                full = cmd + pths + flags
                if pths != "":
                    subprocess.call(full,shell=True)
                pths = ""
                filecnt = 0
        full = cmd + pths + flags
        if pths != "":
            subprocess.call(full,shell=True)
        pths = ""

    return datadir
