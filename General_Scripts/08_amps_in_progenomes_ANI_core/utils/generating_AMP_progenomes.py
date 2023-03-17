import numpy as np
import pandas as pd


def header_maker(val: int) -> str:
    val = str(val).zfill(12)
    val = '_'.join([val[idx:idx+3] for idx in range(0, 12, 3)])
    return f'GMSC10.SMORF.{val}'


def getsample(x, df):
    if x in df.index:
        return df.loc[x, 'sample']
    else:
        return np.nan


df = pd.read_table('../../data_folder/GMSC10.ProGenomes2.coords.txt.gz',
                   sep='.',
                   names=['taxid', 'sample', 'contigloc'])

df['genes'] = [header_maker(x) for x in df.index]

x = df[['genes', 'sample']].set_index('genes')

df = pd.read_table('../../data_folder/GMSC10.Macrel_05.AMPs.tsv.gz')

df['sample'] = df.Access.apply(lambda w: getsample(w, x))

df = df.dropna()

df = df[['sample', 'Sequence']].drop_duplicates()

df.to_csv('data/ref_progenomes_AMPs.tsv.gz',
          sep='\t',
          header=True,
          index=None)
          
