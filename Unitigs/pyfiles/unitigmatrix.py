import pandas as pd
import numpy as np

def get_unitigs(datadir, dname):
    unitigdict = {}
    unitigs = []
    numsamples = 0
    with open(datadir + '/processed_data/' + str(dname) + '_graph.gfa_colored.gfa','r') as f:
        for l in f:
            temp = str.split(l)
            if temp[0] == 'S':
                colors = str.split(temp[4],':')[2]
                numsamples=len(colors)
                colors = np.array(map(int, list(colors)))
                unitigdict[temp[2]] = colors
                unitigs.append(int(temp[1])-1)
    indexes = [x for x in range(numsamples)]

    df = pd.DataFrame(unitigdict,index=indexes)
    df.to_pickle(datadir + '/processed_data/' + str(dname) + '_unitigs.pkl')
    return df
