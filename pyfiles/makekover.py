import os
import pandas as pd
import numpy as np
from collections import Counter

colnames = ['id', 'assembly', 'genus', 'species', 'seqfile', 'cntfile', 'meta']
loadpth = '/home/liam/compare' + '/processed_data/clean.csv'
readcsv = pd.read_csv(loadpth, names=colnames)
paths = readcsv.seqfile.tolist()
inds = readcsv.id.tolist()
species = readcsv.species.tolist()
specdict = {}
index = 0
for s in set(species):
    specdict[s] = index
    index += 1
iddict = {}
meliencoded = []
aborencoded = []
suisencoded = []
for x in range(len(inds)):
    inds[x] = 'genome_' + str(inds[x])
    paths[x] = '/home/liam/compare' + paths[x]
    if species[x] == 'melitensis':
        meliencoded.append(1)
        aborencoded.append(0)
        suisencoded.append(0)
    elif species[x] == 'abortus':
        aborencoded.append(1)
        meliencoded.append(0)
        suisencoded.append(0)
    elif species[x] == 'suis':
        suisencoded.append(1)
        meliencoded.append(0)
        aborencoded.append(0)
    else:
        meliencoded.append(0)
        suisencoded.append(0)
        aborencoded.append(0)

iddict['ids'] = inds
iddict['genomes'] = paths
melimap = {}
melimap['ids'] = inds
melimap['species'] = meliencoded
abormap = {}
abormap['ids'] = inds
abormap['species'] = aborencoded
suismap = {}
suismap['ids'] = inds
suismap['species'] = suisencoded
genomes = pd.DataFrame(iddict)
melimap = pd.DataFrame(melimap)
abormap = pd.DataFrame(abormap)
suismap = pd.DataFrame(suismap)
with open('/home/liam/kovertest/genomes.tsv', 'w') as w:
    w.write(genomes.to_csv(index=False,sep='\t',header=False))
with open('/home/liam/kovertest/melispecies.tsv', 'w') as w:
    w.write(melimap.to_csv(index=False,sep='\t',header=False))
with open('/home/liam/kovertest/aborspecies.tsv', 'w') as w:
    w.write(abormap.to_csv(index=False,sep='\t',header=False))
with open('/home/liam/kovertest/suisspecies.tsv', 'w') as w:
    w.write(suismap.to_csv(index=False,sep='\t',header=False))
    
splits = np.load('/home/liam/compare/processed_data/foldsplits.npy', allow_pickle=True)
for x in range(5):
    with open('/home/liam/kovertest/fold' + str(x+1) + 'train.txt', 'w') as w:
        for i in splits[0][x]:
            w.write(inds[i] + '\n')
        for i in splits[1][x]:
            w.write(inds[i] + '\n')
    with open('/home/liam/kovertest/fold' + str(x+1) + 'test.txt', 'w') as w:
        for i in splits[2][x]:
            w.write(inds[i] + '\n')
            