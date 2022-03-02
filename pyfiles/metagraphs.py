import pandas as pd
import numpy as np
import torch
import os
import sys
from torch_geometric.data import Data

def build_folds():
    datadir = '/home/liam/compare'
    filepth = sys.argv[1]
    outnames = []
    for name in sys.argv[2:]:
        outnames.append(name)
    ####
    # Building the final feature graphs for each individual sample. This is done by building a vector of
    # presence/absence in each sample contained in the pangenome for each unitig (graph node) used to query
    # resulting in a graph of shape [#_nodes][#_pangenome_samples] which is ~655 samples for the current dataset
    ####
    for x in range(5):
        colours = []
        with open(datadir + '/processed_data/fold' + str(x+1) + 'graphsamples.txt', 'r') as f:
            for l in f:
                temp = str.rsplit(l)[0]
                colours.append(temp)
        donegraphs = []
        for outname in outnames:
            ####
            # Building initial graph by building vectors of to/from unitig links contained in graph
            ####
            graphname = datadir + '/processed_data/' + outname + '.gfa'
            with open(graphname, 'r') as f:
                numnodes = 0
                fromarr = []
                toarr = []
                for l in f:
                    temp = str.split(l)
                    if temp[0] == 'S':
                        numnodes = int(temp[1])
                    elif temp[0] == 'L':
                        fromnum = int(temp[1]) - 1
                        tonum = int(temp[3]) - 1
                        if fromnum < 0:
                            fromnum = 0
                        if tonum < 0:
                            tonum = 0
                        fromarr.append(fromnum)
                        toarr.append(tonum)
            g = np.zeros((numnodes,656),dtype=np.float32)
            ####
            # Collecting presence/absence vectors from the queries done in last step
            ####
            graphname = datadir + '/processed_data/fold' + str(x+1) + '_' + outname + '.fasta.search'
            with open(graphname, 'r') as f:
                for l in f:
                        temp = str.split(l)
                        g[int(temp[0])][colours.index(temp[1])] = 1
            fromarr = np.array(fromarr)
            toarr = np.array(toarr)
            edge_index = np.zeros((2,len(toarr)), dtype=np.long)
            edge_index[0] = fromarr
            edge_index[1] = toarr
            y = np.zeros(1,dtype=np.float32)
            y[0] = 5
            y = torch.tensor(y, dtype=torch.long)
            feat = torch.tensor(g, dtype=torch.float32)
            edge_index = torch.tensor(edge_index, dtype=torch.long)
            ####
            # Building Pytorch-geo data object from collected data and adding to list of graphs for this fold
            ####
            data = Data(x=feat, edge_index=edge_index, y=y)
            data.graphind = 1
            donegraphs.append(data)
        print(donegraphs)
        torch.save(donegraphs, datadir + '/processed_data/' + outname + 'fold' + str(x+1) + 'dataset.pkl')
        
    return datadir

build_folds()