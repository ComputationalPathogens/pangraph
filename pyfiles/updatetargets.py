import pandas as pd

def update_targets(datadir, targpth):
    colnames = ['id', 'assembly', 'genus', 'species', 'seqfile', 'cntfile', 'meta']
    data = pd.read_csv(datadir + '/processed_data/counts.csv', names=colnames)
    newtargets = []
    if targpth == '':
        raise SystemExit('Custom targets were selected but no file was provided (use --customtargetspath {path/to/target}')
    with open(targpth,'r') as f:
        for l in f:
            temp = str.rsplit(l)
            newtargets.append(temp)
            
    if len(data['species']) != len(newtargets):
        raise SystemExit('Length of custom targets not equal to number of samples')
    
    data['species'] = newtargets
    data.to_csv(datadir + '/processed_data/counts.csv', index=False, header=False)