import csv
import numpy as np
import pandas as pd
import os
import re
from collections import Counter


def build_index(datadir, dname):

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
    filepth = datadir + '/processed_data/' + str(dname) + '_metadata.csv'
    datapth = datadir + '/samples/' + str(dname) + '/'
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
        ext = ''
        pth, assembly, organism, genus, species, cnt = 'empty', 'empty', 'empty', 'empty', 'empty', 'empty'
        for file in files:
            ext = os.path.splitext(file)[-1].lower()
            if ext == ".fna":
                temp = os.path.basename(subdir)
                name = os.path.splitext(file)[0]
                pth = '/samples/' + str(dname) + '/' + temp + '/' + name + '.fna'
                id += 1
            if ext == ".fasta":
                temp = os.path.basename(subdir)
                name = os.path.splitext(file)[0]
                pth = '/samples/' + str(dname) + '/' + temp + '/' + name + '.fasta'
                id += 1
            if ext == ".gz":
                temp = os.path.basename(subdir)
                extr = os.path.splitext(file)[0]
                name = os.path.splitext(extr)[0]
                fastapth = datapth + temp + '/' +  extr
                pth = '/samples/' + str(dname) + '/' + temp + '/' + name + '.fna'
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