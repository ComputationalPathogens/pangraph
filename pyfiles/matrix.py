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

    genome_row = np.zeros((num_cols), dtype=np.dtype('uint32'))
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

    return genome_row

def build_matrix(datadir, filename = '/processed_data/cleanwcounts.csv'):
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
    if os.path.isfile(datadir + '/processed_data/features.pkl'):
        return datadir
    files_path = datadir + filename
    i = 0
    files = get_file_names(files_path)
    for seq in itertools.product(chars, repeat=11):
        dna = "".join(seq)
        rev = Seq.reverse_complement(dna)
        if dna > rev:
            dna = rev
        if not dna in cols:
            cols[dna] = i
            i += 1
            
    x = np.asarray(files)
    numcols = i
    numrows = len(x)
    kmer_matrix = np.zeros((numrows,numcols),dtype=np.dtype('uint32'))
    rowindex = 0
    with ProcessPoolExecutor(max_workers=None) as ppe:
        for row in ppe.map(get_kmer_counts, files, itertools.repeat(numcols), itertools.repeat(cols), itertools.repeat(datadir)):
            rows[rowindex] = rowindex
            kmer_matrix[rowindex,:] = row
            rowindex += 1
    
    matrixdf = pd.DataFrame(kmer_matrix, columns=cols.keys())
    saves = datadir + '/processed_data/features.pkl'
    matrixdf.to_pickle(saves)
    return datadir


