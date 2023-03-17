import os
import numpy as np
import pandas as pd
import pickle as pkl

from glob import glob
from workfams import families
from openfasta import openfasta
from prot_download import batchdownload
from prot_download import get_assemblies


def checkfile(f):
    '''
    check the integrity of the compressed file f
    '''
    import subprocess
    s = subprocess.run(["gzip",
                        "-t",
                        f],
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE)
    return s.returncode
         
         
print('Load Progenomes2 clusters')
clusters = pd.read_table('data/proGenomes2.1_specI_clustering.tab',
                         sep='\t',
                         header=None,
                         names=['specI', 'sample'])

clusters['sample'] = clusters['sample'].apply(lambda x: x.split('.')[-1])
clusters = clusters.groupby('specI').apply(lambda x: x['sample'].tolist())
clusters = clusters.reset_index().rename({0: 'samples'}, axis=1)
clusters['L'] = clusters.samples.apply(lambda x: len(x))
clusters = clusters[clusters.L >= 10].sort_values(by='L').drop('L', axis=1)

# check for already done clusters
alreadydone = []
for f in glob('analysis/*_proteins.tsv.gz'):
    f = f.split('/')[-1]
    f = '_'.join(f.split('_')[0:3])
    alreadydone.append(f)

if os.path.exists('Notes.txt'):
    with open('Notes.txt', 'rt') as f:
        for row in f:
            row = row.replace('ERROR for cluster ', '')
            row = row.split(',')[0]
            alreadydone.append(row)

# clean to start the process
for f in glob('*.gz'): os.remove(f)
    
with open('Notes.txt', 'a') as output:
    for _, c, i in clusters.itertuples():
        if c in alreadydone:
            print(f'Already analyzed {c}, passing to next')
        else:
            print(f'Working with cluster {c} -------------')
            batchdownload(i)
            for f in glob('*.gz'):
                if checkfile(f):
                    x = f.split('.')[0]
                    links = get_assemblies(x,
                                   download=True)
                    if checkfile(f):
                        print('.')
                    else:
                        os.remove(f)
            N = len(list(glob('*.gz'))) 
            if N >= 10:
                seqs = dict()
                for f in glob('*.gz'):
                    hs = openfasta(f)
                    for h, s in hs:
                        if s in seqs:
                            seqs[s].append(h)
                        else:
                            seqs[s] = [h]
                    os.remove(f)
                nlen, L = [], len(seqs)
                for k, v in seqs.items():
                    nseqs = [x.split('_')[0] for x in v]
                    nseqs = set(nseqs)
                    nseqs = len(nseqs)
                    nlen.append(nseqs*100/N)
                nlen = np.asarray(nlen)
                res = []
                for cutoff in range(0, 101):
                    k = sum(nlen >= cutoff)
                    res.append((cutoff, k, k*100/L))
                df = pd.DataFrame(res,
                                  columns=['cutoff_prevalence',
                                           'number of proteins',
                                           'percent of unique proteins'])
                df.to_csv(f'analysis/{c}_proteins.tsv.gz',
                          sep='\t',
                          header=True,
                          index=None)
                with open(f'analysis/{c}_seqs.pkl', 'wb') as handle:
                    pkl.dump(seqs,
                             handle,
                             protocol=pkl.HIGHEST_PROTOCOL)    
                families()
            else:
                print(f'ERROR for cluster {c}, no sequences appended')
                output.write(f'ERROR for cluster {c}, {N} download genomes, no sequences appended\n')
            for f in glob('*.gz'):
                os.remove(f)
                             
