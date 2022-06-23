import csv
import numpy as np
import pandas as pd
import os
import re
from collections import Counter


def build_metadata(datadir, filename = '/processed_data/metadata.csv'):

    """
    Parameters
    ----------
    datadir : Base directory of nextflow execution

    Returns
    -------
    datadir : Base directory of nextflow execution
    
    Disc
    ----
    Generates a .csv to record the Assembly#, Genus, Species and file paths
    for complete genome sequences/kmer counts

    """
    filepth = datadir + filename
    datapth = datadir + '/refseq/bacteria/'
    f = open(filepth, 'w', newline='')
    writer = csv.writer(f)
    id = -1
    
    meta = re.compile("Date|Submitter|Assembly method|Genome coverage|Sequencing technology")
    #change over to use pandas write_csv
    for subdir, dirs, files in os.walk(datapth):
        if re.search("mers", str(subdir)):
            continue
        dirs.sort()
        files.sort()
        metadict = {}
        pth, assembly, organism, genus, species, cnt = 'empty', 'empty', 'empty', 'empty', 'empty', 'empty'
        for file in files:
            ext = os.path.splitext(file)[-1].lower()
            if ext == ".fna" or ext == '.fasta':
                temp = os.path.basename(subdir)
                name = os.path.splitext(file)[0]
                if ext == '.fna':
                        pth = '/refseq/bacteria/' + temp + '/' + name + '.fna'
                else:
                        pth = '/refseq/bacteria/' + temp + '/' + name + '.fasta'
            if ext == ".gz":
                temp = os.path.basename(subdir)
                extr = os.path.splitext(file)[0]
                name = os.path.splitext(extr)[0]
                fastapth = datapth + temp + '/' +  extr
                pth = '/refseq/bacteria/' + temp + '/' + name + '.fna'
                cmd = 'gunzip ' + fastapth
                os.system(cmd)
            if ext == ".txt":
                fp = subdir + '/' + file
                with open(fp, 'r') as f:
                    lines = f.readlines()
                assembly = lines[0]
                assembly = assembly[18:]
                organism = lines[1]
                organism = organism[18:]
                species = organism.split(" ")
                genus = species[0]
                species = species[1]
                if species == 'sp.':
                    species = "unknown"
                for l in lines[2:]:
                    if meta.search(l):
                        metadict[str(meta.search(l)[0])] = str(re.split(': +', str(meta.split(l)[1]))[1])[:-1]
                        
        if id > -1:
            writer.writerow([id, assembly, genus, species, pth, cnt, metadict])
        id += 1

    return datadir

def clean_outliers(k, datadir, filename = '/processed_data/metadata.csv'):
    """
    Parameters
    ----------
    k : The amount of folds for k-fold validation, will remove any classes with
    less than k occurences so stratified k-fold can be used
    
    datadir : Base directory of nextflow execution

    Returns
    -------
    datadir : Base directory of nextflow execution
    
    Disc
    ----
    'Cleans' data by removing any samples with less than k occurences so data
    can be used in stratified k-fold

    """
    k = int(k)
    colnames = ['id', 'assembly', 'genus', 'species', 'seqfile', 'cntfile', 'meta']
    csvpth = datadir + filename
    data2 = pd.read_csv(csvpth, names=colnames)
    species = data2.species.tolist()
    labels = np.asarray(species)
    count = Counter(labels)
    brucunknown = []

    for index, row in data2.iterrows():
        if row['species'] == 'unknown' and row['genus'] == 'Brucella':
            brucunknown.append(row['seqfile'])
        if count[row['species']] < 10 or row['species'] == 'unknown':
            data2.drop(index, axis=0, inplace=True)
    data2.reset_index(drop=True,inplace=True)
    newinds = []
    with open(datadir + '/processed_data/unknownset.txt', 'w') as f:
        for file in brucunknown:
            f.write(datadir + file + '\n')
    for x in range(len(data2.index)):
        newinds.append(x)
    newinds = {'id':newinds}
    newinds = pd.DataFrame(newinds)
    data2['id'] = newinds['id']
    cleanpth = datadir + '/processed_data/clean.csv'
    data2.to_csv(cleanpth, index=False, header=False)
    return datadir
