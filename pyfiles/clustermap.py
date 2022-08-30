from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.preprocessing import LabelEncoder
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import sys

def create(datadir, dname):
    sys.setrecursionlimit(100000)
    colnames=['id','assembly','genus','species','seqfile','cntfile', 'meta']
    data = pd.read_pickle(datadir+'/processed_data/' + str(dname) + '_featuresfiltered.pkl')
    labels = pd.read_csv(datadir+'/processed_data/' + str(dname) + '_counts.csv', names=colnames)
    labels = labels.species.tolist()
    data.index = [x for x in range(len(data))]
    label_enc = LabelEncoder()
    label_enc = label_enc.fit(labels)
    clrs = sns.color_palette("husl",14)
    lut = dict(zip(label_enc.classes_, clrs))
    rowclrs = []
    for x in labels:
        rowclrs.append(lut[x])
    """
    sel = SelectKBest(f_classif,k=50000)
    kbestdata = sel.fit_transform(data,labels)
    kbestmask = sel.get_support()
    kbestlabels = data.columns[kbestmask]
    kbestdf = pd.DataFrame(kbestdata, columns=kbestlabels)
    kbestdf.index = [x for x in range(len(data))]
    """
    kbestdf = data.sample(n=50000,axis='columns')
    snsplot = sns.clustermap(kbestdf, cmap='flare',row_colors=rowclrs, vmin = 0, vmax = 5, method='complete',row_cluster = True,col_cluster=False)
    print(snsplot)
    handles = [Patch(facecolor=lut[name]) for name in lut]
    snsplot.savefig(datadir+"/processed_data/" + str(dname) + "_ClusterMap.png")
    return snsplot
    plt.legend(handles, lut, title='Species',
               bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right')
    
create('/home/liam/Phylogeny','flaviv')
    

