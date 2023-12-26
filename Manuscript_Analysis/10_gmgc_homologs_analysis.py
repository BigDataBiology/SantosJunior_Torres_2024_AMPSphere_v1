#!/usr/bin/env python

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from collections import Counter
from itertools import combinations
from scipy.stats import kruskal, mannwhitneyu

plt.rcParams['svg.fonttype'] = 'none'

M8_NAMES = ['query', 'target', 'evalue',
           'gapopen', 'pident', 'nident',
           'qstart', 'qend', 'qlen',
           'tstart', 'tend', 'tlen',
           'alnlen', 'raw', 'bits',
           'cigar', 'qseq', 'tseq',
           'qheader', 'theader', 'qaln',
           'taln', 'qframe', 'tframe',
           'mismatch', 'qcov', 'tcov']

gmgc = pd.read_table('../data_folder/result_gmgc.m8.xz',
                     names=M8_NAMES,
                     header=None,
                     usecols=['query', 'evalue', 'bits', 'pident', 'qstart', 'qend', 'qlen', 'tstart', 'tend', 'tlen'])

# creating some measures
# make zero-based
gmgc['tstart'] -= 1
# fix for macrel initial methionine deletion
gmgc.loc[gmgc.tstart == 1, 'tstart'] = 0
gmgc['pct_start'] = gmgc.eval('tstart * 100 / tlen')
gmgc['pct_end']   = gmgc.eval('tend * 100 / tlen')
gmgc['pct_amp']   = gmgc.eval('(1 + qend - qstart)*100 / qlen')


# histogram of pct_start
gmgc.sort_values(
    by=['evalue', 'bits', 'pident'],
    ascending=[True, False, False],
    inplace=True,
)

df = gmgc.groupby('query').head(1)

# df.to_csv('filtered_gmgc_homologs.tsv.gz', sep='\t', header=True, index=None)

fig,ax = plt.subplots()
sns.histplot(data=df,
            x='pct_start',
            bins=100,
            color='black',
            ax=ax)

ax.set_xlabel('Match start (% of target length)')
ax.set_ylabel('Counts')

f = df.eval('pct_start <= 25 | pct_start >= 75').mean()
print(f'{f:.1%} of hits begin in the initial or final 25% of target protein')


# ## Enrichment of ortholog groups among the homologs

data = pd.read_table('../data_folder/adjust_significant_function.csv.xz',
                     sep='\t',
                     header='infer')

data = data[['eggnog_OG',
             'count_AMP',
             'count_GMGC',
             'times',
             'p_adjust']]

data.columns = ['EggNOG ortholog group',
                'Counts in the homologs of c_AMPs',
                'Counts in the redundant GMGC',
                'Enrichment (fold)',
                'Adjusted P-value']

e5 = pd.read_table('../data_folder/e5_annotations.tsv.gz',
                  sep='\t',
                  header=None)

e5.columns = ['tax',
              'EggNOG ortholog group',
              'class',
              'function']

e5 = e5[['EggNOG ortholog group',
         'class',
         'function']].drop_duplicates()

e5.drop('function', axis=1, inplace=True)
e5 = e5.drop_duplicates()
k = e5['EggNOG ortholog group'].value_counts()
k = k[k==1]
e5 = e5[e5['EggNOG ortholog group'].isin(k.index)]

print('Working with OGs')

x = data.merge(on='EggNOG ortholog group',
               right=e5)

cmap = {'S': ('Unknown function', (43/255, 142/255, 112/255, 255/255)),
        'L': ('Replication, recombination,\nand repair', (190/255, 99/255, 29/255, 255/255)),
        'J': ('Translation, ribosomal\nstructure, and biogenesis', (124/255, 120/255, 171/255, 255/255)),
        'K': ('Transcription', (101/255, 149/255, 47/255, 255/255)),
        'C': ('Energy production\nand conversion', (207/255, 65/255, 138/255, 255/255)),
        '-': ('Other', (202/255, 157/255, 30/255, 255/255))}

x['class'] = x['class'].apply(lambda x: x if x in cmap else '-')
x['class'] = x['class'].apply(lambda x: cmap.get(x)[0])

a = x['class'].value_counts() * 100 / len(x)
print(a)

a, b = min(x['Enrichment (fold)']), max(x['Enrichment (fold)'])
print(f'Minimum enrichment was {a}-fold, and the maximum was {b}-fold')

cutoff = 1000
n = len(x[x['Enrichment (fold)'] < cutoff])
n = n*100/len(x)
print(f'{n:.2f}% of OGs were enriched less than {cutoff}-fold')

x.rename({'class': 'Functional class'}, axis=1, inplace=True)

x = x[x['Enrichment (fold)'] >= 1.5]
x = x[x['Adjusted P-value'] < 0.05]

x.to_csv('outputs/suppTableS3.tsv',
         sep='\t',
         header=True,
         index=None)


a = x.groupby('Functional class').agg('size')
b = x[x['Counts in the homologs of c_AMPs'] <= 10].groupby('Functional class').agg('size')
c = x[x['Counts in the homologs of c_AMPs'] <= 5].groupby('Functional class').agg('size')
d = pd.concat([a, b, c], axis=1)
d.columns = ['c > 10', '5 < c <= 10', 'c <= 5']
d *= -1
d = d.reset_index()
d = d.melt(id_vars='Functional class')
order = ['Unknown function',
         'Translation, ribosomal\nstructure, and biogenesis',
         'Replication, recombination,\nand repair',
         'Other',
         'Energy production\nand conversion',
         'Transcription']
d = pd.concat([d[d['Functional class'] == w] for w in order])


colors = {'c > 10': (140/255, 45/255, 4/255),
          '5 < c <= 10': (254/255, 153/255, 41/255),
          'c <= 5': (255/255, 247/255, 188/255)}

for w in colors:
    sns.barplot(data=d[d['variable'] == w],
                y='Functional class',
                x='value',
                color=colors[w],
                label=w)

ax.legend()
fig.tight_layout()
fig.savefig('outputs/GMGC_og_functions.svg')


x = x.query('`Counts in the homologs of c_AMPs` >= 20')

order = x.groupby('Functional class')
order = order['Enrichment (fold)'].quantile([0.5])
order = order.sort_values(ascending=False).index
order = [y[0] for y in order]

cmap = {'Unknown function': (43/255, 142/255, 112/255, 255/255),
        'Translation, ribosomal\nstructure, and biogenesis': (190/255, 99/255, 29/255, 255/255),
        'Replication, recombination,\nand repair': (124/255, 120/255, 171/255, 255/255),
        'Other': (202/255, 157/255, 30/255, 255/255),
        'Energy production\nand conversion': (207/255, 65/255, 138/255, 255/255),
        'Transcription': (101/255, 149/255, 47/255, 255/255)}

fig, ax = plt.subplots()
sns.boxplot(data=x,
            y='Functional class',
            x='Enrichment (fold)',
            color='white',
            width=0.4,
            showfliers=False,
            order=order,
            ax=ax)

sns.stripplot(data=x,
              y='Functional class',
              x='Enrichment (fold)',
              s=3,
              palette=cmap,
              order=order,
              ax=ax)

ax.set_xscale('log')
ax.legend('')
fig.tight_layout()
fig.savefig('outputs/gmgc_functional_enrichment.svg')

# ## Statistical testing of enrichment and abundance

print('Test of c_AMP homolog counts for each OG')
c = x.groupby('Functional class').apply(lambda w: w['Counts in the homologs of c_AMPs'].tolist())
_, p = kruskal(*c)
print(f'Kruskal result: p={p:.2e}')

print('Paired:')
print('category1\tcategory2\tP-value')
for i, j in combinations(c.index, 2):
    _, pp = mannwhitneyu(c[i], c[j])
    if pp < 0.05: print(f'{i}\t{j}\t{pp}')

print('Test of c_AMP enrichment')
c = x.groupby('Functional class').apply(lambda w: w['Enrichment (fold)'].tolist())
_, p = kruskal(*c)
print(f'Kruskal result: p={p:.2e}')

print('Paired:')
print('category1\tcategory2\tP-value')
for i, j in combinations(c.index, 2):
    _, pp = mannwhitneyu(c[i], c[j])
    if pp < 0.05: print(f'{i}\t{j}\t{pp}')

