#!/usr/bin/env python
# coding: utf-8

import pickle
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from collections import Counter
from scipy.stats import pearsonr

from environments import color_map
color_map = {k:tuple(list(v)+[.5]) for k, v in color_map.items()}

def calc_rare(infile):
    print('# counting nsamples')
    df = pd.read_table(infile)
    m = df['nsamples'].mean()
    n = sum(df['nsamples'] < m)
    q50 = df['nsamples'].median()
    print(f'The median number of samples in which a AMP is detected is {q50}')
    print(f'The mean number of samples in which a AMP is detected is {m}')
    print(f'{n} AMPs were detected less than or equal to the mean number of detections across all AMPSphere candidates')
    print(f'This represents {n*100/len(df):.2f}% of all AMPs')


def plot_test_rare(infile, envdict, ofile):
    cutoffs = [0.01, 0.1, 1, 10]
    envdict['total'] = sum(envdict.values())
    df = pd.read_table(infile, sep='\t', header='infer')
    df = df.set_index('AMP')
    df['total'] = df.sum(axis=1)
    headers = ['Habitat',
               'Pearson_R',
               'P-value',
               'N samples',
               'Total AMPs detected',
               '0.01%_samples',
               '0.1%_samples',
               '1%_samples',
               '10%_samples']
    with open(f'outputs/{ofile}.tsv', 'w') as handle:
        handle.write('\t'.join(headers)+'\n')
        for col in df.columns:
            if envdict[col] >= 100:
                ncs = [int(envdict[col]*c/100) for c in cutoffs]
                totamps = len(df[df[col] > 0])
                a = df.loc[df[col] > 0, col].value_counts()
                a = a.sort_index()
                a = a.reset_index()
                a.columns = ['samples', 'amps']
                a['1/k'] = 1/a['samples']
                r, p = pearsonr(a.amps, a['1/k'])
                n = []
                for c in ncs:
                    i = a.loc[a['samples'] <= c, 'amps'].sum()
                    n.append(i)
                n = [str(x) for x in n]
                vals = [col, f'{r:.3f}', f'{p:.2E}',
                        str(envdict[col]),
                        str(totamps),
                        n[0], n[1], n[2], n[3]]
                handle.write('\t'.join(vals)+'\n')
                if col in color_map:
                    sns.scatterplot(y=a['amps'],
                                    x=a['samples'],
                                    label=col,
                                    alpha=0.25,
                                    color=color_map.get(col),
                                    s=3)
                else:
                    sns.scatterplot(y=a['amps'],
                                    x=a['samples'],
                                    label=col,
                                    alpha=0.25,
                                    s=3)
        plt.yscale('log')
        plt.xscale('log')
        plt.xlabel('Detections\n(Number of samples)')
        plt.ylabel('c_AMPs')
        plt.tight_layout()
        plt.savefig(f'outputs/{ofile}.svg')
        plt.close()
        _ = envdict.pop('total')



with open('../data_folder/envo.pkl', 'rb') as handle:
    envo = pickle.load(handle)
h = dict(Counter(envo['general_envo_name'].values()))
hi = dict(Counter(envo['high'].values()))


calc_rare('../data_folder/unique_map_cAMPs_nsamples.tsv.xz')

plot_test_rare('../data_folder/unique_map_cAMPs_high_habs.tsv.xz',
               hi,
               'output_hi')

plot_test_rare('../data_folder/unique_map_cAMPs_general_habs.tsv.xz',
               h,
               'output_gen')

