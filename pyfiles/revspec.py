import ahocorasick
import pandas as pd
from Bio import Seq, SeqIO, SearchIO, SeqRecord
from datetime import datetime

datadir = '/home/liam/compare'

colnames = ['id', 'assembly', 'genus', 'species', 'seqfile', 'cntfile', 'meta']
samples = pd.read_csv(datadir + '/processed_data/clean.csv', names=colnames)
seqfile = samples.seqfile.tolist()
species = samples.species.tolist()
print(samples.shape)
for file, spec in zip(seqfile,species):
    if spec == 'melitensis':
        with open(datadir + file, 'r') as f:
            seqs = []
            for record in SeqIO.parse(f, 'fasta'):
                rseq = record.seq[::-1]
                seqs.append(rseq)
            newrow = [str(len(samples)),str(len(samples)) + '_rev', 'Brucella', 'rev_melitensis', file+'_REV.fna','empty','empty']
            samples.loc[len(samples)] = newrow
            seqreqs = []
            idcount = 0
            for s in seqs:    
                seqreq = SeqRecord.SeqRecord(seq=s,id=str(len(samples)-1)+'_rev_' + str(idcount))
                idcount += 1
                seqreqs.append(seqreq)
            with open(datadir + file + '_REV.fna', 'w') as output:
                SeqIO.write(seqreqs, output, 'fasta')
        
samples.to_csv('/home/liam/compare/revtest.csv')
        