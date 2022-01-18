import ahocorasick
import pandas as pd
from Bio import Seq, SeqIO
from datetime import datetime

datadir = '/home/liam/compare'
primers = ['GACGAACGGAATTTTTCCAATCCC', 'AAATCGCGTCCTTGCTGGTCTGA', 'CGGGTTCTGGCACCATCGTCG', 'GCGCGGTTTTCTGAAGGTTCAGG', 'TGCCGATCACTTAAGGGCCTTCAT']
labels = ['BA', 'BM', 'BO', 'BS', 'REV']
lengths = {'BA':500,'BM':731,'BO':976,'BS':285}
AMOS = ahocorasick.Automaton()
for p,l in zip(primers,labels):
    AMOS.add_word(p,l)
    rev = Seq.reverse_complement(p)
    lab = 'RC_' + l
    AMOS.add_word(rev,lab)
    
AMOS.make_automaton()

colnames = ['id', 'assembly', 'genus', 'species', 'seqfile', 'cntfile', 'meta']
samples = pd.read_csv(datadir + '/processed_data/clean.csv', names=colnames)
seqfile = samples.seqfile.tolist()
species = samples.species.tolist()

seqs = []
for file in seqfile:
    with open(datadir + file, 'r') as f:
        seq = ''
        for l in f:
            temp = str.rsplit(l)[0]
            if temp[0] != '>':
                seq += temp
        seqs.append(seq)

print(datetime.now())
for seq, specs in zip(seqs,species):
    last_index, last_value = 0, ''
    first = True
    val_hits = []
    for end_index, original_value in AMOS.iter_long(seq):
        start_index = end_index - len("GACGAACGGAATTTTTCCAATCCC") + 1
        if not first:
            left = str.split(original_value,'_')
            right = str.split(last_value, '_')
            if right[0] == 'RC' and left[0] != 'RC':
                dist = start_index - last_index
                if right[-1] != left[-1]: #not both REV primers
                    if right[-1] != 'REV': #get length of target gene based on primer that isnt REV
                        if lengths[right[-1]] + 100 > dist:
                            val_hits.append((dist,left,right))
                    if left[-1] != 'REV':
                        if lengths[left[-1]] + 100 > dist:
                            val_hits.append((dist,left,right))
            elif right[0] != 'RC' and left[0] == 'RC':
                dist = start_index - last_index
                if right[-1] != left[-1]: #not both REV primers
                    if right[-1] != 'REV': #get length of target gene based on primer that isnt REV
                        if lengths[right[-1]] + 100 > dist:
                            val_hits.append((dist,left,right))
                    if left[-1] != 'REV':
                        if lengths[left[-1]] + 100 > dist:
                            val_hits.append((dist,left,right))
        
        last_index = start_index
        last_value = original_value
        first = False
    print("SPECIES: " + specs)
    print("valid hits:")
    for hit in val_hits:
        print(hit)
print(datetime.now())