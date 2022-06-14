import itertools
import pandas as pd
import os
from Bio import Seq, SeqIO
import numpy as np
from concurrent.futures import ProcessPoolExecutor

def get_file_names(filepth):
    """
    Parameters
    ----------
    filepth : Path to clean data .csv to pull kmer count file paths

    Returns
    -------
    files : List of file paths for counted kmers

    """

    colnames = ['id', 'assembly', 'genus', 'species', 'seqfile', 'cntfile', 'meta']
    data = pd.read_csv(filepth, names=colnames)
    files = data.cntfile.tolist()
    return files


def get_kmer_counts(filename, num_cols, col_index, datadir):
    """
    Parameters
    ----------
    filename : Path to file of counted kmers
    num_cols : Amount of features
    col_index : Dictionary mapping of Kmer:Count for building feature row

    Returns
    -------
    genome_row : Completed feature row of kmer counts for input sample

    """
    #count each kmer that exists, look into jellyfish doing this for us, collect the highest freq kmers?
    genome_row = np.zeros((num_cols), dtype=np.dtype('uint32'))
    presence_row = np.zeros((num_cols), dtype=np.dtype('uint32'))
    append = datadir + filename
    #Same method used in Computational-pathogens/acheron
    with open(append) as f:
        for record in SeqIO.parse(f, 'fasta'):
            kmercount = record.id
            kmercount = int(kmercount)
            seq = record.seq
            seq = str(seq)
            if len(seq) < 11:
                print("Invalid key: %s" % (seq))
                print(filename)
            index = col_index[seq]
            genome_row[index] = kmercount
            if kmercount > 0:
                presence_row[index] = 1
            else:
                presence_row[index] = 0

    return genome_row, presence_row

def build_matrix(datadir, filename = '/processed_data/counts.csv'):
    """
    Parameters
    ----------
    datadir : Base directory of nextflow execution

    Returns
    -------
    datadir : Base directory of nextflow execution
    
    Desc
    ----
    Builds feature matrix of [NumSamples]*[NumFeatures] size and saves it 
    to features.pkl
    
    """
    cols = {}
    rows = {}
    chars = "ACGT"

    files_path = datadir + filename
    i = 0
    if not os.path.isfile(datadir + '/processed_data/features.pkl'):
        files = get_file_names(files_path)
        isnum = 0
        for f in files:
            with open(datadir + f, 'r') as counts:
                isnum = 0
                for l in counts:
                    if isnum % 2 != 0:
                        seq = str.rsplit(l)[0]
                        rev = Seq.reverse_complement(seq)
                        if seq > rev:
                            seq = rev
                        if seq not in cols:
                            cols[seq] = i
                            i += 1
                    isnum += 1
            """         < 31mer code
            for seq in itertools.product(chars, repeat=31):
                dna = "".join(seq)
                rev = Seq.reverse_complement(dna)
                if dna > rev:
                    dna = rev
                if not dna in cols:
                    cols[dna] = i
                    i += 1
            """

        
        
    if not os.path.isfile(datadir + '/processed_data/features.pkl'):
        x = np.asarray(files)
        numcols = i
        numrows = len(x)
        kmer_matrix = np.zeros((numrows,numcols),dtype=np.dtype('uint32'))
        pres_matrix = np.zeros((numrows,numcols),dtype=np.dtype('uint32'))
        rowindex = 0
        with ProcessPoolExecutor(max_workers=None) as ppe:
            for row, pres in ppe.map(get_kmer_counts, files, itertools.repeat(numcols), itertools.repeat(cols), itertools.repeat(datadir)):
                rows[rowindex] = rowindex
                kmer_matrix[rowindex,:] = row
                pres_matrix[rowindex,:] = pres
                rowindex += 1
        
        matrixdf = pd.DataFrame(kmer_matrix, columns=cols.keys())
        saves = datadir + '/processed_data/features.pkl'
        matrixdf.to_pickle(saves)
        presdf = pd.DataFrame(pres_matrix, columns=cols.keys())
        pressaves = datadir + '/processed_data/featurespresence.pkl'
        presdf.to_pickle(pressaves)
        
        colsums = presdf.sum(axis=0)
        filtered = []
        for s in colsums:
            if colsums[s] >= 5:
                filtered.append(s)
        filtersaves = datadir + '/processed_data/featuresfiltered.pkl'
        filtereddf = presdf.filter(filtered, axis=1)
        filtereddf.to_pickle(filtereddf)
    return datadir


