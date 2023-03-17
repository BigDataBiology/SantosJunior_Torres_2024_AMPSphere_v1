import os
import numpy as np
import pandas as pd
import pickle as pkl

from glob import glob
from subprocess import run
from nuc_download import batchdownload
from nuc_download import get_assemblies
from res_ani_test import checkres

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

print('Load AMPSphere AMPs')

amps = pd.read_table('data/complete_gmsc_pgenomes_metag.tsv.gz', sep='\t', header='infer')
amps = amps[amps.is_metagenomic == False]
amps['sample'] = amps['sample'].apply(lambda x: '.'.join(x.split(".")[1:]))
amps = amps[['amp', 'sample']].drop_duplicates()

print('Load AMPSphere families')

spheres = pd.read_table('data/SPHERE_v.2022-03.levels_assessment.tsv.gz',
                        sep='\t', header='infer')

spheres = spheres[['AMP accession', 'SPHERE_fam level III']]
spheres.columns = ['amp', 'fam']
spheres = spheres.set_index('amp').to_dict()
spheres = spheres['fam']

amps['fam'] = amps['amp'].apply(lambda x: spheres.get(x, 'NA'))

print('Checking not to do')
alreadynot = []
if os.path.exists('Notes.txt'):
    with open('Notes.txt', 'rt') as f:
        for row in f:
            row = row.replace('ERROR for cluster ', '')
            row = row.split(',')[0]
            alreadynot.append(row)

print('Clean to start the process')
for f in glob('*.gz'): os.remove(f)
    
print('Start your engines')
for _, c, i in clusters.itertuples():
    print(f'Working with cluster {c}')
    if os.path.exists(f'analysis/{c}_ANI.tsv'):
        alreadynot.append(c)
        
    if c in alreadynot:
        print(f'Skipping {c}, passing to next')
    else:
        batchdownload(i)
        for f in glob('*.fna.gz'):
            if checkfile(f):
                x = f.split('.')[0]
                links = get_assemblies(x,
                               download=True)
                if checkfile(f):
                    print('.')
                else:
                    os.remove(f)

        N = len(list(glob('*.fna.gz')))  
        
        with open('temporary_list.txt', 'wt') as of: 
            for x in list(glob('*.fna.gz')): of.write(f'{x}\n')                    
        
        if (N >= 10):          
            run(['fastANI',
                 '--ql', 'temporary_list.txt',
                 '--rl', 'temporary_list.txt',
                 '-o', f'analysis/{c}_ANI.tsv',
                 '-t', '4'])
            os.remove('temporary_list.txt')
#            checkres(c, amps)
              
    print('Clean to start the process')
    for f in glob('*.fna.gz'): os.remove(f)
                        
