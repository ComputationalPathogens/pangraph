from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.preprocessing import LabelEncoder
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import sys

def create(datadir):
    sys.setrecursionlimit(100000)
    colnames=['id','assembly','genus','species','seqfile','cntfile', 'meta']
    data = pd.read_pickle(datadir+'/processed_data/featuresfiltered.pkl')
    labels = pd.read_csv(datadir+'/processed_data/counts.csv', names=colnames)
    labels = labels.species.tolist()
    data.index = labels
    label_enc = LabelEncoder()
    label_enc = label_enc.fit(labels)
    clrs = sns.color_palette("husl",14)
    lut = dict(zip(label_enc.classes_, clrs))
    rowclrs = []
    for x in labels:
        rowclrs.append(lut[x])
    
    sel = SelectKBest(f_classif,k=10000)
    kbestdata = sel.fit_transform(data,labels)
    kbestmask = sel.get_support()
    kbestlabels = data.columns[kbestmask]
    kbestdf = pd.DataFrame(kbestdata, columns=kbestlabels)
    kbestdf.index = labels
    snsplot = sns.clustermap(kbestdf, row_colors=rowclrs, vmin = 0, vmax = 10, method='complete')
    handles = [Patch(facecolor=lut[name]) for name in lut]
    plt.legend(handles, lut, title='Species',
               bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right')
    
    snsplot.savefig(datadir+"/processed_data/ClusterMap.png")
    return snsplot

