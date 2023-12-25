#!/usr/bin/env python

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
from random import shuffle
from scipy.stats import norm, shapiro
from itertools import chain, permutations

from environments import higher_level, color_map, animal_guts

data = pd.read_table('data/gmsc_amp_genes_envohr_source.tsv.gz',
                     sep='\t',
                     header='infer')

# filter duplicates
data = data.query('is_metagenomic')
data = data[['amp', 'general_envo_name']].drop_duplicates()
data = data.groupby('general_envo_name')['amp'].apply(set)

# select environments with at least 100 peptides
data = data[data.map(len) >= 100]


# add the high_level environment
data = data.reset_index()
data['high'] = [higher_level.get(x, 'other') for x in data.general_envo_name]
data.set_index('general_envo_name', inplace=True)

# ### For low level habitats
# calculate overlap
df = []
combs = permutations(data.index, 2)
for i, j in combs:
    set_i = data.loc[i, 'amp']
    set_j = data.loc[j, 'amp']
    n = len(set_i.intersection(set_j))
    df.append((i, j, n))

# including doubled pair
for i in data.index:
    df.append((i,
               i,
               len(data.loc[i, 'amp'])))

# formatting table
df = pd.DataFrame(df, columns=['env1', 'env2', 'overlap'])
df = df.pivot(index='env1', columns='env2', values='overlap')

# normalize
df = df * 100 / df.max(axis=0)

# organize environments
df['high'] = [higher_level.get(x, 'other') for x in df.index]
# sort high level habitat and low level envo names
df = df.sort_index().sort_values(by='high')
df = df.drop('high', axis=1)
# sort columns by high level habitat
ncol = df.index
df = df[ncol]


# create a color map for the general envo having
# as keys the high level environment
colors = [color_map[higher_level[i]] for i in df.index]

g = sns.clustermap(data=df,
                   col_colors=colors,
                   row_colors=colors,
                   cmap='YlOrBr',
                   row_cluster=False,
                   col_cluster=False)

for label in color_map:
    g.ax_col_dendrogram.bar(0,
                            0,
                            color=color_map[label],
                            label=label,
                            linewidth=0)

g.ax_col_dendrogram.legend(loc="center", ncol=6)
g.cax.set_position([.15, .2, .03, .45])


# ### For high level environments

df = data.groupby('high')['amp'].apply(lambda x: set(chain.from_iterable(x)))


newdf = []
for i, j in permutations(df.index, 2):
    set_i = df.loc[i]
    set_j = df.loc[j]
    n = len(set_i.intersection(set_j))
    newdf.append((i, j, n))

for i in df.index:
    newdf.append((i,
                  i,
                  len(df.loc[i]))
                )

newdf = pd.DataFrame(newdf, columns=['env1', 'env2', 'overlap'])
newdf = newdf.pivot(index='env1', columns='env2', values='overlap')


# normalize
newdf = newdf * 100 / newdf.max(axis=0)

# create mask of zeros
mask = np.zeros_like(newdf)
mask[np.tril_indices_from(mask)] = True


fig, ax = plt.subplots()
sns.heatmap(newdf.astype('int'),
            annot=False,
            cmap="YlOrBr",
            mask=mask, square=True,
            ax=ax)
fig.tight_layout()


# ### For animal guts
# getting AMP overlap from guts in each host


newdf = []
for i, j in permutations(animal_guts, 2):
    set_i = data.loc[i, 'amp']
    set_j = data.loc[j, 'amp']
    n = len(set_i.intersection(set_j))
    newdf.append((i, j, n))

for i in animal_guts:
    newdf.append((i,
                  i,
                  len(data.loc[i, 'amp']))
                )

newdf = pd.DataFrame(newdf, columns=['env1', 'env2', 'overlap'])
newdf = newdf.pivot(index='env1', columns='env2', values='overlap')


# normalize
newdf = newdf * 100 / newdf.max(axis=0)

# create mask of zeros
mask = np.zeros_like(newdf)
mask[np.tril_indices_from(mask)] = True

fig, ax = plt.subplots()
sns.heatmap(newdf.astype('int'),
            annot=False,
            cmap="YlOrBr",
            mask=mask, square=True,
            ax=ax)

fig.tight_layout()


# ### For human body sites

human_body = ['human skin', 'human respiratory tract',
              'human mouth',
              'human digestive tract', 'human gut',
              'human urogenital tract']

# we need to account for human mouth and human saliva
# they will become just human mouth
x = data.loc['human saliva', 'amp'].union(data.loc['human mouth', 'amp'])
data = data.drop(['human mouth', 'human saliva'], axis=0)
data.loc['human mouth'] = [x, 'other human']


# calculating overlaps
newdf = []
for i, j in permutations(human_body, 2):
    set_i = data.loc[i, 'amp']
    set_j = data.loc[j, 'amp']
    n = len(set_i.intersection(set_j))
    newdf.append((i, j, n))

# adding pairs of same index
for i in human_body:
    newdf.append((i,
                  i,
                  len(data.loc[i, 'amp']))
                )

# formatting result
newdf = pd.DataFrame(newdf, columns=['env1', 'env2', 'overlap'])
newdf = newdf.pivot(index='env1', columns='env2', values='overlap')


# normalize
newdf = newdf * 100 / newdf.max(axis=0)

# create mask of zeros
mask = np.zeros_like(newdf)
mask[np.tril_indices_from(mask)] = True

fig,ax = plt.subplots()
sns.heatmap(data=newdf.astype('int'),
            annot=False,
            cmap="YlOrBr",
            mask=mask, square=True,
            ax=ax)
fig.tight_layout()


