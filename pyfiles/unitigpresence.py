import pandas as pd
import numpy as np
from Bio import Seq, SeqIO



def get_presence(datadir):
    #splits = np.load(datadir + '/processed_data/foldsplits.npy', allow_pickle=True)
    #graph_loader = DataLoader(graphs, batch_size=1)
    colnames = ['id', 'assembly', 'genus', 'species', 'seqfile', 'cntfile', 'meta']
    samples = pd.read_csv(datadir + '/processed_data/clean.csv', names=colnames)
    species = samples.species.tolist()


    graphids = samples.id.tolist()


    currind = 0

    baseunitigs = {}

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
                        baseunitigs[seq] = currind
                        currind += 1
                lnum += 1
    tig_matrix = np.zeros((len(graphids),currind),dtype=np.dtype('uint32'))
    for i in graphids:
        genome_row = np.zeros((currind), dtype=np.dtype('uint32'))
        with open(datadir + '/processed_data/fasta/graph' + str(i) + '.fasta', 'r') as file:
            lnum = 0
            for l in file:
                if lnum % 2 != 0:
                    seq = str.rsplit(l)[0]
                    rev = Seq.reverse_complement(seq)
                    if seq > rev:
                        seq = rev
                    index = baseunitigs[seq]
                    genome_row[index] = 1
                lnum += 1
            tig_matrix[i,:] = genome_row
    matrixdf = pd.DataFrame(tig_matrix,columns=baseunitigs.keys())
    print(matrixdf)