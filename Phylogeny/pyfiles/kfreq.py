import numpy as np
import pandas as pd
import subprocess
import os
from concurrent.futures import ProcessPoolExecutor
import itertools

def count_kmer(file, rootdir, ksize = 11):
    """
    Parameters
    ----------
    rootdir : Base directory of nextflow execution

    Returns
    -------
    merpth : the path for the generated mer_counts.jf file that will be used for stats later

    """
    merpth = ''
    id = 0
    ext = os.path.splitext(str(file))
    dirname = os.path.dirname(file)
    if ext[-1] == ".fna":
        try:
            mkpth = rootdir + dirname + '/' + str(ksize) + '-freqmers/'
            os.mkdir(mkpth)
        except OSError:
            pass
            
        merpth = rootdir + dirname + '/' + str(ksize) + '-freqmers/' + 'mer_counts.jf'
        seqpth = rootdir + file
        cmd = 'jellyfish count -m ' + str(ksize) + ' -s 100M -t 10 -C -o ' + merpth + ' ' + seqpth
        os.system(cmd)
        id = id + 1
    else:
        print("error, no .fna file exists at given directory")
        id = id + 1
        
    countout = kmer_freq(merpth, ksize)
    
    return countout

def kmer_freq(merpth, ksize):

    args = ['jellyfish','stats',merpth]
    output = subprocess.check_output(args,text=True)
    output = output.split("\n")
    stats = []
    for s in output:
        stats.append(s.split(" ")[-1])
    freq = int(stats[2])/int(stats[1])
    print(stats[3])

    return freq

def genome_freqs(datadir, ksize):
    
    colnames = ['id', 'assembly', 'genus', 'species', 'seqfile', 'cntfile', 'meta']
    filepth = datadir + '/processed_data/clean.csv'
    data = pd.read_csv(filepth, names=colnames)
    files = data.seqfile.tolist()

    genome_matrix = np.zeros(len(files))
    genomeindex = 0
    with ProcessPoolExecutor(max_workers=None) as ppe:
        for countjf in ppe.map(count_kmer, files, itertools.repeat(datadir), itertools.repeat(ksize)):
            genome_matrix[genomeindex] = countjf
            genomeindex += 1
            
    outpth = datadir + '/processed_data/' + str(ksize) + '-merfreq.npy'
    np.save(outpth, genome_matrix)
    return outpth