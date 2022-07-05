import matplotlib.pyplot as plt
import numpy as np
import math

def freqBoxPlot(datadir):
    labels = []
    log2freqs = []
    for x in range(11,32,2):
        labels.append(str(x) + '-mer')
        temp = []
        log2 = []
        temp = np.load(str(datadir) + '/processed_data/' + str(x) + '-merfreq.npy')
        for y in temp:
            log2.append(math.log2(y))
        log2freqs.append(log2)
        
    plt.boxplot(log2freqs, vert=True, labels=labels, whis=99999)
    plt.suptitle('K-mer Frequency Among Genome Dataset', fontsize=12)
    plt.ylabel('Log2(Frequency)')
    plt.xlabel('Kmer Length')
    plt.xticks(fontsize=6)
    plt.yscale('log',base=2)
    plt.savefig("KmerFreq.png")
