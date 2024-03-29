import pandas as pd
import sys
import ast
import scipy.cluster
from sklearn.feature_selection import SelectKBest, f_classif

def build_data(datadir,dname):
    data = pd.read_pickle(datadir + '/processed_data/' + str(dname) + '_featuresfiltered.pkl')
    colnames = ['id','assembly','genus','species','seqfile','cntfile','meta']
    metadata = pd.read_csv(datadir+'/processed_data/' + str(dname) + '_counts.csv', names=colnames)
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
    sel = SelectKBest(f_classif,k=1000000)
    kbestdata = sel.fit_transform(data,labels)
    kbestmask = sel.get_support()
    kbestlabels = data.columns[kbestmask]
    kbestdf = pd.DataFrame(kbestdata, columns=kbestlabels, dtype=float)
    kbestdf.index = labels
    return kbestdf, meta

def hierarchal(data, k):
    
    link = scipy.cluster.hierarchy.linkage(data, method = 'complete')
    clusters = scipy.cluster.hierarchy.cut_tree(link, n_clusters=k)
    
    return clusters

def kmeans(data, k):
    centroids, clusters = scipy.cluster.vq.kmeans2(data, k=k, minit='points')
    
    return clusters

def comp_clusters(datadir, numclust,dname):
    sys.setrecursionlimit(100000)
    data,meta = build_data(datadir,dname)
    print("Starting kClustering")
    kClust = kmeans(data,numclust)    
    with open(datadir + '/processed_data/' + str(dname) + '_kclust.txt', 'w') as f:        
        for h in kClust:
            f.write(str(h) + '\n')
    print("Starting hClustering")
    hClust = hierarchal(data, numclust) 
    with open(datadir + '/processed_data/' + str(dname) + '_hclust.txt', 'w') as f:        
        for h in hClust:
            f.write(str(h) + '\n')

    

comp_clusters('/home/liam/Phylogeny', 6, 'anthrax')
