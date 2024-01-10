#!/usr/bin/env python
# coding: utf-8
import random
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from tqdm import tqdm
from scipy.stats import norm
from scipy.stats import shapiro
from scipy.stats import pearsonr, spearmanr

import environments

data = pd.read_table('../data_folder/gmsc_amp_genes_envohr_source.tsv.gz')

# eliminate genes/amps without environment (from ProGenomes)
data = data[~data.general_envo_name.isna()]

# attribute high-level habitat to genes/amps
data['high'] = data.general_envo_name.map(lambda x: environments.higher_level.get(x, 'other'))


# testing multihabitat AMPs
def get_multihabitat(df, level=None):
    '''
    counts the number of multi-habitat AMPs (present in at least 3 environments)

    :inputs:
    data frame containing at least AMP, general_envo_name, high-level habitat
    level which stats for the habitat or high-level environment

    :outputs:
    length of the list of multi-habitat c_AMPs
    '''
    if level == None:
        level = 'high'
    if level == 'low':
        level = 'general_envo_name'

    h = df[['amp', level]].drop_duplicates()
    h = h.amp.value_counts()
    h = h[h > 1].index
    return h


l = get_multihabitat(data, 'high')
l0 = get_multihabitat(data, 'low')
print(f'Multi-habitat AMPs: {len(l0):,}\nMulti-high-level-habitat: {len(l):,}')


k = data[['amp', 'high']].drop_duplicates()['amp'].value_counts()
k = k[k>1].index

with open('outputs/multihabitat_highlevel.txt', 'wt') as out:
    for i in k: out.write(f'{i}\n')

k = data[['amp', 'general_envo_name']].drop_duplicates()['amp'].value_counts()
k = k[k>1].index

with open('outputs/multihabitat_generalenvo.txt', 'wt') as out:
    for i in k: out.write(f'{i}\n')


# ### Permutation test
#
# We shuffle the high-level habitat annotation for the samples, and then, calculate the number of multi-habitat c_AMPs. This operation is repeated 100 times, and then we calculate the average and standard deviation of the distribution of random results. Using Shapiro-Wilk test, we check if the random distribution is normal, and if it is, we calculate the Z-score for the result obtained for AMPSphere. The Z-score is then converted into a p-value to support our conclusions.

# testing significance
def shuffle_test(df, level=None):
    if level == None: level = 'high'
    elif level == 'low': level = 'general_envo_name'

    habitat = df.set_index('sample')[level].to_dict()

    values = [v for _, v in habitat.items()]
    random.shuffle(values)

    for k, v in zip(habitat, values):
        habitat[k] = v

    altdf = df.copy()
    altdf[level] = altdf['sample'].map(lambda x: habitat.get(x))

    return altdf


#test high level habitats
test = []
for _ in tqdm(range(100)):
    altdf = shuffle_test(data, 'high')
    altdf = len(get_multihabitat(altdf, 'high'))
    test.append(altdf)

print(test)

_, p = shapiro(test)
print(f'The Shapiro-Wilks test returned p = {p}')

avg, std = np.mean(test), np.std(test)

print(f'Average number of random multi-habitat AMPs - {avg}, with std = {std}')

Z = (len(l) - avg) / std

pz = norm.sf(abs(Z))

print(f'The number of multi-habitat AMPs in AMPSphere was {len(l)}')
print(f'It was {Z} * std of the random distribution')
print(f'This gives us a p-value of {pz}')


#test habitats
test_low = []
for _ in tqdm(range(100)):
    altdf = shuffle_test(data, 'low')
    altdf = len(get_multihabitat(altdf, 'low'))
    test_low.append(altdf)

print(test_low)

_, p = shapiro(test_low)

print(f'The Shapiro-Wilks test returned a p = {p}')

avg, std = np.mean(test_low), np.std(test_low)

print(f'Average number of random multi-habitat AMPs - {avg}, with std = {std}')

Z = (len(l0) - avg) / std

pz = norm.sf(abs(Z))

print(f'The number of multi-habitat AMPs in AMPSphere was {len(l0):,}')
print(f'It was {Z} * std of the random distribution')
print(f'This gives us a p-value of {pz}')

