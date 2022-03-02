import pandas as pd
import numpy as np
import torch
import torch.nn.functional as F
from Bio import Seq, SeqIO
from torch.nn import Linear
from torch_geometric.nn import (GraphConv, SAGPooling, global_mean_pool, JumpingKnowledge)
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import precision_recall_fscore_support, accuracy_score
from captum.attr import Saliency, IntegratedGradients


device = torch.device('cuda')

def model_forward(node_mask, data, model):
    integrated = Data(x=node_mask, edge_index=data.edge_index, y=data.y, batch=data.batch, ptr = data.ptr, graphind = data.graphind)
    out = model(integrated)
    return out

def explain(data, model, target=0, exp_type='ig'):
    input_mask = torch.ones((data.x.shape[0],data.x.shape[1])).to(device)
    input_mask = torch.tensor(input_mask)
    if exp_type == 'ig':
        ig = IntegratedGradients(model_forward)
        mask = ig.attribute(input_mask, target=data.y[0],additional_forward_args=(data,model,),internal_batch_size=data.x.shape[0])
    elif exp_type == 'sal':
        sal = Saliency(model_forward)
        mask = sal.attribute(input_mask, target=data.y[0],additional_forward_args=(data,model,))
    node_mask = np.abs(mask.cpu().detach().numpy())
    if node_mask.max() > 0:
        node_mask = node_mask / node_mask.max()
    return node_mask

def importance(datadir, graphs, model, split, curspecies):
    features = []
    splits = np.load(datadir + '/processed_data/foldsplits.npy', allow_pickle=True)
    graph_loader = DataLoader(graphs, batch_size=1)
    colnames = ['id', 'assembly', 'genus', 'species', 'seqfile', 'cntfile', 'meta']
    samples = pd.read_csv(datadir + '/processed_data/clean.csv', names=colnames)
    species = samples.species.tolist()
    enc = 0
    indspec = {}
    specdict = {}
    labels_unencoded = []
    #TURN INTO HELPER FUNCTION IN DIFFERENT CLASS EVERYWHERE THIS EXISTS
    for s in set(species):
        labels_unencoded.append(s)
    labels_unencoded.sort()
    for s in labels_unencoded:
        specdict[s] = 0
    graphids = samples.id.tolist()
    speciesmask = []
    specimportance = {}
    for l in labels_unencoded:
        specimportance[l] = np.zeros(657, dtype=int)
    currind = 0
    for graph in splits[1][split]:
        specimportance[species[graph]][currind] = 1
        currind += 1
        """
    for graph in splits[1][split]:
        if species[graph] == curspecies:
            speciesmask.append(1)
        else:
            speciesmask.append(0)
        """
    print(specimportance)
    #speciesmask.append(0)
    baseunitigs = {}
    indunitigs = {}
    unitiginds = {}
    curind = 0
    for i in graphids:
        with open(datadir + '/processed_data/fasta/graph' + str(i) + '.fasta', 'r') as file:
            lnum = 0
            for l in file:
                if lnum % 2 != 0:
                    seq = str.rsplit(l)[0]
                    rev = Seq.reverse_complement(seq)
                    if seq > rev:
                        seq = rev
                    if seq not in baseunitigs:
                        baseunitigs[seq] = 0
                        indunitigs[curind] = seq
                        unitiginds[seq] = curind
                        
            
                        curind += 1
                lnum += 1
    
    pair_importance = {}
    for d in graph_loader:
        unitigs = baseunitigs.copy()
        d = d.to(device)
        expl = explain(d, model, 'sal')
        importance = np.zeros(expl.shape[0])
        maximp = 0
        maxind = 0
        maxclass = ""
        for node in range(len(expl)):
            impdict = specdict.copy()
            for l in labels_unencoded:
                impdict[l] = expl[specimportance[l]].sum()
            #for hit in range(len(expl[node])):
            #    if speciesmask[hit]:
            #        importance[node] += expl[node][hit]
            #if importance[node] > maximp:
            #    maximp = importance[node]
            #    maxind = node
            for l in labels_unencoded:
                if impdict[l] > maximp:
                    maximp = impdict[l]
                    maxclass = l
        
        print(maximp)
        print(maxclass)
        imprankings = sorted(impdict.items(), key=lambda x:x[1], reverse=True)
        print(imprankings)
        count = 0
        seqs = []
        with open(datadir + '/processed_data/fasta/graph' + str(d.graphind.item()) + '.fasta') as f:
            for l in f:
                if count % 2 != 0:
                    seq = str.rsplit(l)[0]
                    rev = Seq.reverse_complement(seq)
                    if seq > rev:
                        seq = rev
                    seqs.append(seq)
                count += 1
        for s in range(len(seqs)):
            unitigs[seqs[s]] += importance[s]
        with open(datadir + '/processed_data/graphs/graph' + str(d.graphind.item()) + '.gfa','r') as f:
            for l in f:
                temp = str.rsplit(l)
                if temp[0] == 'L':
                    seqL = int(temp[1]) - 1
                    seqR = int(temp[3]) - 1
                    pair = (unitiginds[seqs[seqL]], unitiginds[seqs[seqR]])
                    avgimp = ((importance[seqL] + importance[seqR]) / 2)
                    if pair not in pair_importance:
                        pair_importance[pair] = avgimp
                    else:
                        pair_importance[pair] += avgimp


    #pairsort = dict(sorted(pair_importance.items(), key=lambda item: item[1]), reverse=True)
    pairlist = sorted(pair_importance.items(), key=lambda x:x[1], reverse=True)
    print("#############NEIGHBOR SEQUENCES###############")
    print(pairlist[0])
    print(pairlist[0][0], pairlist[0][1])
    print(indunitigs[int(pairlist[0][0][0])], indunitigs[int(pairlist[0][0][1])])
    
    print(pairlist[1])
    print(pairlist[1][0], pairlist[1][1])
    print(indunitigs[int(pairlist[1][0][0])], indunitigs[int(pairlist[1][0][1])])
    
    maximp = 0
    maxseq = ''
    imp = []
    for k,v in unitigs.items():
        imp.append((v,k))
    imp.sort(reverse = True)
    #print(imp[0:5])
    for k in unitigs.keys():
        if unitigs[k] > maximp:
            maximp = unitigs[k]
            maxseq = k
    print("########SINGLE SEQUENCES#########")
    print(maxseq)
    print(maximp)
    return maxseq