def cluster_by_ani(infile: str):
    import os
    import pickle as pkl
    import gzip
    import pandas as pd
    
    from random import sample
    
    cluster = '_'.join(infile.split('_')[:-1])
    
    data = pd.read_table(infile,
                         sep='\t',
                         header=None,
                         names=['s1', 's2', 'ANI', 'm', 't'])
    
    data['s1'] = data.s1.apply(lambda x: x.replace('.fna.gz', ''))
    data['s2'] = data.s2.apply(lambda x: x.replace('.fna.gz', ''))
    
    x = data[data.ANI > 99.5]
    x = x.groupby('s1')
    x = x.apply(lambda w: set(w.s2.tolist()))
    x = x.reset_index()
    x = x.rename({'s1': 'provk',
                  0: 'provc'}, axis=1)
    listofsets = []
    for k in x.provc:
        if k in listofsets: pass
        else: listofsets.append(k)
    
    representatives = dict()
    for k in listofsets:
        kk = sample([*k], k=1)
        for element in k:
            representatives[element] = kk[0]
    
    with open(f'strains_{cluster}.pkl', 'wb') as handle:
        pkl.dump(representatives,
                 handle,
                 protocol=pkl.HIGHEST_PROTOCOL)    
    
    with open(f'strains_{cluster}.txt', 'w') as ofile:
        for w in set(representatives.values()):
            ofile.write(w+'\n')   
    
