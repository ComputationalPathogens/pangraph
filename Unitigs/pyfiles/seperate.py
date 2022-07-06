import os

def seperate_viral(datadir, filepth, dname):
    try:
        os.mkdir(datadir + '/samples/' + str(dname) + '/')
    except:
        pass
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
                    mkpth = datadir + '/samples/' + str(dname) + '/' + str(currsample)
                    os.mkdir(mkpth)
                except:
                    pass
                w = open(datadir + '/samples/' + str(dname) + '/' + str(currsample) + '/' + str(currsample) + '.fasta', 'w')
                w.write(l)
            else:
                w.write(l)
