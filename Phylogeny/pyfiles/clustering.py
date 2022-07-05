import pandas as pd
import sys
import ast
import scipy.cluster
from sklearn.feature_selection import SelectKBest, f_classif
from collections import Counter
from sklearn.feature_selection import VarianceThreshold

def build_data(datadir):
    data = pd.read_pickle(datadir + '/processed_data/featuresfiltered.pkl')
    colnames = ['id','assembly','genus','species','seqfile','cntfile','meta']
    metadata = pd.read_csv(datadir+'/processed_data/counts.csv', names=colnames)
    labels = metadata.species.tolist()
    """
    sk_obj = VarianceThreshold(threshold=0.75)
    fitdata = sk_obj.fit_transform(data)
    print(fitdata.shape)
    meta = []
    metalist = metadata.meta.tolist()
    for x in metalist:
        meta.append(ast.literal_eval(x))
    fitdata = fitdata.astype(float)
    """
    meta = []
    metalist = metadata.meta.tolist()
    for x in metalist:
        meta.append(ast.literal_eval(x))
    sel = SelectKBest(f_classif,k=100000)
    kbestdata = sel.fit_transform(data,labels)
    kbestmask = sel.get_support()
    kbestlabels = data.columns[kbestmask]
    kbestdf = pd.DataFrame(kbestdata, columns=kbestlabels)
    kbestdf.index = labels
    return kbestdf, meta

def hierarchal(data, k):
    
    link = scipy.cluster.hierarchy.linkage(data, method = 'complete')
    clusters = scipy.cluster.hierarchy.cut_tree(link, n_clusters=k)
    
    return clusters

def kmeans(data, k):
    centroids, clusters = scipy.cluster.vq.kmeans2(data, k=k, minit='points')
    
    return clusters

def comp_clusters(datadir, numclust):
    sys.setrecursionlimit(100000)
    data,meta = build_data(datadir)
    
    hClust = hierarchal(data, numclust)
    
    index = 0
    

    for h in hClust:
        print("HClust#: ", h)
        index += 1
    
    with open(datadir + '/processed_data/hclust2.txt', 'w') as f:        
        for h in hClust:
            f.write('Hierarchal Cluster #' + str(h) + '\n')
            print("Hierarchal Cluster #", h)
    
    
    

comp_clusters('/home/liam/Phylogeny', 25)
