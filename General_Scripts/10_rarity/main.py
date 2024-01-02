import os
import lzma
import pickle
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import cm

from tqdm import tqdm
from glob import glob
from os.path import exists
from itertools import groupby
from collections import Counter
from scipy.stats import shapiro, norm, pearsonr
from scipy.sparse import coo_matrix, save_npz, load_npz


def convname(ampcode: str) -> str:
    '''
    Convert AMP access code into row index
    Example: AMP10.000_000 = 0
    '''
    return int(ampcode[6:])


def preprocess_files():
    '''
    Join mapping results for AMP genes against sample reads
    Outputs:
    samples_list = list of samples which the index in the list match to the column index in the sparse matrix
    data = sparse matrix (lil) of presence/absence with the rows being AMPs and columns being samples
    '''
    print('Load dict of genes to proteins')
    headers = pd.read_table('data/headers.tsv.xz',
                            sep='\t',
                            header='infer',
                            index_col='gene')
    n = len(glob('data/uniques/*.xz'))
    tamps = len(set(headers.AMP))
    headers['AMP'] = headers.AMP.apply(convname)
    headers = headers.to_dict()
    headers = headers['AMP']
    print('Processing mapping results')
    samples_list = []
    row, col = [], []
    for idx, data in tqdm(enumerate(glob('data/uniques/*.xz')), total=n):
        data = pd.read_table(data,
                             comment='#',
                             index_col=0)
        samples_list.append(data.columns[0])
        if '-1' in data.index: data.drop('-1', axis=0, inplace=True)
        elif -1 in data.index: data.drop(-1, axis=0, inplace=True)
        if len(data) != 0:
            amps = {headers.get(x) for x in data.index}
            for x in amps:
                row.append(x)
                col.append(idx)
    print('Producing sparse matrix')
    print(f'Matrix def: ({n}, {tamps})')
    data = np.ones_like(row, dtype=np.uint8)
    data = coo_matrix((data, (row, col)),
                      shape=(tamps, n),
                      dtype=np.uint8)
    save_npz('analysis/sparse_matrix.npz', data)  # to load: scipy.sparse.load_npz('sparse_matrix.npz')
    with open('analysis/saved_sample_cols.txt', 'w') as handle:
        for x in samples_list: handle.write(f'{x}\n')
    return (samples_list, data)


def load_precomputed(override: bool=False):
    if not override:
        print('Trying to retrieve saved files info')
        if exists('analysis/saved_sample_cols.txt'):
            print('Loading samples')
            samples = []
            for row in open('analysis/saved_sample_cols.txt', 'r'):
                samples.append(row.strip())
            if exists('analysis/sparse_matrix.npz'):
                print('Loading sparse_matrix')
                df = load_npz('analysis/sparse_matrix.npz')
                return (samples, df)
    print('Starting general process')
    return preprocess_files()


def metadata():
    '''
    Import metadata
    Outputs:
    1 = dictionary with keys as samples and values as the general environments
    2 = dictionary with keys as samples and values as the high-level environments
    '''
    with open('data/envo.pkl', 'rb') as handle:
        envo = pickle.load(handle)
    return (envo['general_envo_name'], envo['high'])


def export_tables(protmap_res: list, samples_h: list, samples_g: list):
    from collections import Counter
    print('Creating and setting variables')
    nsamples, ghabs, hhabs = [], [], []
    g, h = set(samples_g), set(samples_h)
    a, b = np.array(samples_h), np.array(samples_g)
    print('Creating tables for output')
    ghabs_out = lzma.open('analysis/unique_map_cAMPs_general_habs.tsv.xz', 'wt')
    gnames = '\t'.join(['AMP', *tuple(g)])
    ghabs_out.write(f'{gnames}\n')
    nsamples_out = lzma.open('analysis/unique_map_cAMPs_nsamples.tsv.xz', 'wt')
    nsamples_out.write('AMP\tnsamples\n')
    hhabs_out = lzma.open('analysis/unique_map_cAMPs_high_habs.tsv.xz', 'wt')
    hnames = '\t'.join(['AMP', *tuple(h)])
    hhabs_out.write(f'{hnames}\n')
    print('Processing results')
    for idx in tqdm(range(protmap_res.shape[0]), total=protmap_res.shape[0]):
        v = protmap_res[idx].nonzero()[1]
        n = str(len(v))
        name = str(idx).zfill(6)
        name = 'AMP10.'+'_'.join([name[0:3], name[3:]])
        nsamples_out.write(f'{name}\t{n}\n')
        g_coo = dict(Counter(b[v]))
        g_coo = [ g_coo.get(x, 0) for x in g ]
        g_coo.insert(0, name)
        g_coo = '\t'.join([str(w) for w in g_coo])
        ghabs_out.write(f'{g_coo}\n')
        h_coo =  dict(Counter(a[v]))
        h_coo = [ h_coo.get(x, 0) for x in h ]
        h_coo.insert(0, name)
        h_coo = '\t'.join([str(w) for w in h_coo])
        hhabs_out.write(f'{h_coo}\n')
    print('Closing inputs')
    hhabs_out.close()
    nsamples_out.close()
    ghabs_out.close()


samples, df = load_precomputed()
h, hi = metadata()
samples_g = [h.get(x) for x in samples]
samples_h = [hi.get(x) for x in samples]
df = df.tocsr()
print('Export tables')
export_tables(df, samples_h, samples_g)

