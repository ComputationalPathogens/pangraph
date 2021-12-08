import pandas as pd
import numpy as np
import os
import subprocess

def check_kmers(datadir):
    filepath = datadir + '/processed_data/clean.csv'
    colnames = ['id', 'assembly', 'genus', 'species', 'seqfile', 'cntfile', 'meta']
    dumpname = 'mer_counts_dumpsbf.fa'
    
    data = pd.read_csv(filepath, names=colnames)
    id = 0
    data['seqfile'] = data['seqfile'].astype('str')
    files = data.seqfile.tolist()
    species = data.species.tolist()
    indspec = {}
    specdict = {}
    labels_unencoded = []
    for s in set(species):
        labels_unencoded.append(s)
    labels_unencoded.sort()
    for s in labels_unencoded:
        specdict[s] = 0
    checkmers = [
    'AGGGGAATAAGGGAATAGGGGAGTATGTTTG',
'CATCGCCGAAAGCGGGCATGAACCGCTGTCC',
'GACGCTTGAACCGAGGGATGCAAATGTTGGA',
'CGTTTCACACTTTTGCTGGAAATGCTCTATA',
'TCCACTGGACTGATTTCTGATCCCGCTTCGA',
'GATCGACGCCAATGCCAAGGGCGTTGCCGAC',
'TCTGCATTCAACGTAACCAGATCATAGCGCA',
'CCATTGGCAAGATCAAGGCGATATTCTTCCG',
'GCGTAGCGGTTTTGCGTTGGATAATGCGACC',
'TTCTGGCCAATGGCAGCCCTATTGTCGGCAA',
'CGGCATTCAGCGCCAGAAGGTTGGTCTGGAA',
'GCCAAGGGCGTTGCCGACAACAAGACTGCCA',
'CGACGCAAAACCGCTTCGCACTTTTGGCTCG',
'CCCGAAAACCGTTTCACACTTTTCGGGATGC',
'ATAGAGCATTTCCAGCAAAAGTGTGAAACGG',
'GGCGGGCCATAGCGTTCGGGAACTTCCGCCC',
'GTCATATGCGGATAAAGCGCATAGGACTGGA',
'AACGTACGAAGCCCCACGTTAACATCGGCAC',
'CAGCAAAAGTGTGAAACGGTTTTGCGTTGGA',
'CAAAACCGCTTCGCACTTTTGCTGGAAATGC']

    """
    with open('/home/liam/kovertest/meliresults/model_rule_1_equiv.fasta', 'r') as f:
        linenum = 0
        for l in f:
            if linenum % 3 == 1:
                seq = str.rsplit(l)[0]
                checkmers.append(seq)
            linenum += 1
    """
    #checktig = 'GCCTGATTTGTTGCTCGAACTCTTTTCCGAAGAAATCCCCGCTCTGCCGTCCGCCCGCAGCGAAGCGAGGACGGACAAATTGAATAATGTGCATGCGGAAGGCGTACGCTGATGCCTGATTTGTTGCTCGAACTCTTTTCC'
    #for x in range(len(checktig) - 30):
    #    checkmers.append(checktig[x:x+31])

    ind = 0
    with open(datadir + '/processed_data/check.fasta', 'w') as f:
        for c in checkmers:
            f.write('>' + str(ind) + '\n')
            f.write(c + '\n')
            ind += 1
    checkdict = {}

    
    res = specdict.copy()
    ind = 0
    for spec, seqfile in zip(species, files):
        print("FILE# " + str(ind))
        ind += 1
        dirname = os.path.dirname(seqfile)
        #cmd = 'jellyfish query ' + datadir + dirname + '/31-mers/mer_countsbf.jf ' + mer
        cmd3 = 'jellyfish dump countresults.jf'
        cmd2 = 'jellyfish count -m 31 -s 100M -C -t 16 -o countresults.jf --if processed_data/check.fasta ' + datadir + seqfile
        output = subprocess.check_output(cmd2, shell=True)
        output = subprocess.check_output(cmd3, shell=True)
        output = output.decode('UTF-8')
        output = str.rsplit(output, '\n')
        res[spec] += 1
        for l in range(len(checkmers)):
            num = str.split(output[l*2], '>')
            seq = output[(l*2) + 1]
            if seq not in checkdict:
                checkdict[seq] = specdict.copy()
            if int(num[1]) > 0:
                checkdict[seq][spec] += 1
            
            #output = str.split(output)[1]
            #if int(output) > 0:
            #    res[spec] += 1
            #print(res)
    print(specdict)
    for k in checkdict.keys():
        print(k)
        print(checkdict[k])