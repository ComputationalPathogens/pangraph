import pandas as pd
import numpy as np
import torch
import os
from torch_geometric.data import Data

def build_folds(datadir, dname):
    
    ####
    # Helper Function?
    ####
    splits = np.load(datadir + '/processed_data/' + str(dname) + '_foldsplits.npy', allow_pickle=True)
    colnames = ['id', 'assembly', 'genus', 'species', 'seqfile', 'cntfile', 'meta']
    samples = pd.read_csv(datadir + '/processed_data/' + str(dname) + '_clean.csv', names=colnames)
    species = samples.species.tolist()
    specdict = {}
    indspec = {}
    labels_unencoded = []
    enc = 0
    for s in set(species):
        labels_unencoded.append(s)
    labels_unencoded.sort()
    for s in labels_unencoded:
        specdict[s] = enc
        indspec[enc] = s
        enc += 1
    
    ####
    # Building the final feature graphs for each individual sample. This is done by building a vector of
    # presence/absence in each sample contained in the pangenome for each unitig (graph node) used to query
    # resulting in a graph of shape [#_nodes][#_pangenome_samples] which is ~655 samples for the brucella dataset
    ####
    for x in range(5):
        colours = []
        newgraphc = 0
        graphfeatlength = 0
        with open(datadir + '/processed_data/' + str(dname) + '_fold' + str(x+1) + 'graphsamples.txt', 'r') as f:
            for l in f:
                temp = str.rsplit(l)[0]
                colours.append(temp)
                graphfeatlength += 1
        if not os.path.isfile(datadir + '/processed_data/' + str(dname) + '_fold' + str(x+1) + 'dataset.pkl'):
            donegraphs = []
            shufgraphs = []
            foldsplits = np.concatenate((splits[0][x],splits[2][x]), axis=0)
            ####
            # Building initial graph by building vectors of to/from unitig links contained in graph
            ####
            for graph in foldsplits:
                graphseqs = []
                graphname = datadir + '/processed_data/' + str(dname) + '_graphs/graph' + str(graph) + '.gfa'
                with open(graphname, 'r') as f:
                    numnodes = 0
                    fromarr = []
                    toarr = []
                    for l in f:
                        temp = str.split(l)
                        if temp[0] == 'S':
                            numnodes = int(temp[1])
                            graphseqs.append(temp[2])
                        elif temp[0] == 'L':
                            fromnum = int(temp[1]) - 1
                            tonum = int(temp[3]) - 1
                            if fromnum < 0:
                                fromnum = 0
                            if tonum < 0:
                                tonum = 0
                            fromarr.append(fromnum)
                            toarr.append(tonum)
                g = np.zeros((numnodes,graphfeatlength),dtype=np.float32)
                ####
                # Collecting presence/absence vectors from the queries done in last step
                ####
                graphname = datadir + '/processed_data/' + str(dname) + '_fold' + str(x+1) + '/querygraph_graph' + str(graph) + '.fasta.search'
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
                y[0] = specdict[species[graph]]
                y = torch.tensor(y, dtype=torch.long)
                feat = torch.tensor(g, dtype=torch.float32)
                edge_index = torch.tensor(edge_index, dtype=torch.long)
                ####
                # Building Pytorch-geo data object from collected data and adding to list of graphs for this fold
                ####
                data = Data(x=feat, edge_index=edge_index, y=y)
                data.graphind = int(graph)
                donegraphs.append(data)
                """
                #Code for shuffled species experiment
                if species[graph] == 'melitensis':
                    alty = np.zeros(1,dtype=np.float32)
                    alty[0] = enc
                    alty = torch.tensor(alty, dtype=torch.long)
                    indices = np.arange(g.shape[0])
                    np.random.shuffle(indices)
                    g = g[indices]
                    graphseqs = np.array(graphseqs)
                    graphseqs = graphseqs[indices]
                    altfeat = torch.tensor(g, dtype=torch.float32)
                    altdata = Data(x=altfeat, edge_index = edge_index, y=alty)
                    altdata.graphind = len(species) + newgraphc
                    newgraphc += 1
                    shufgraphs.append(altdata)
                    with open(datadir + '/processed_data/shuffled/graph' + str(int(altdata.graphind)) + '.fasta', 'w') as newf:
                        lines = 0
                        for seqs in graphseqs:
                            newf.write('>' + str(lines) + '\n')
                            newf.write(seqs + '\n')
                            lines += 1
                """
            torch.save(donegraphs, datadir + '/processed_data/' + str(dname) + '_fold' + str(x+1) + 'dataset.pkl')
            #torch.save(shufgraphs, datadir + '/processed_data/fold' + str(x+1) + 'shufdataset.pkl')
    return datadir
