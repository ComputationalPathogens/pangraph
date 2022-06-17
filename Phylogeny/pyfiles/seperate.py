import os

def seperate_viral(datadir, filepth):
    with open(filepth, 'r') as f:
        currsample = 0
        for l in f:
            if l[0] == '>':
                try:
                    w.close()
                    currsample += 1
                except:
                    pass
                try:
                    mkpth = datadir + '/refseq/virus/' + str(currsample)
                    os.mkdir(mkpth)
                w = open(datadir + '/refseq/virus/' + str(currsample) + '/' + str(currsample) + '.fasta', 'w')
                w.write(l)
            else:
                w.write(l)
                
                