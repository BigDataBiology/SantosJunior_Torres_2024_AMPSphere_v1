#!/usr/bin/env python

# ## Testing significance of c_AMPs overlap in habitats
# 
# We tested overlap of c_AMP contents of habitats presenting AMPs in at least 100 samples. Then, a permutation test was performed by shuffling the labels of samples and recalculating the overlap between each pair of habitats 32 times. The average and standard deviation overlap c_AMPs was calculated and the Z-score of the actual measure was computed. The p-value is then calculated using the survival function of *scipy.stats*.

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
from random import shuffle
from scipy.stats import norm, shapiro
from itertools import chain, combinations

from environments import higher_level, color_map, animal_guts

data = pd.read_table('data/gmsc_amp_genes_envohr_source.tsv.gz',
                     sep='\t',
                     header='infer')
data = data.query('is_metagenomic')
data = data[['amp', 'sample', 'general_envo_name']].drop_duplicates()


# eliminating environments with less than 100 samples
df = data[['sample', 'general_envo_name']].drop_duplicates()
ec = df.general_envo_name.value_counts()
envs100 = ec[ec >= 100].index

df = data[data.general_envo_name.isin(set(envs100))].reset_index(drop=True)
df['high'] = df.general_envo_name.map(higher_level.get)

def permtest(df, n: int):
    '''
    Permutation test
    It accepts a data frame (df) consisting of at least three columns: amp, sample, environment
    [ 'general_envo_name' represents the column name with the environment labels ]
    Then it selects subsets of the dataframe by environments v1 and v2
    Shuffles the sample labels and calculate overlap
    '''
    tests, testsh = [], []
    envo, envo_h = set(df.general_envo_name), set(df.high)
    habitat_list, hhabitat_list = list(combinations(envo, 2)), list(combinations(envo_h, 2))
    print('Performing pair-wise comparisons')
    alist = df.groupby('general_envo_name').apply(lambda x: set(x.amp))
    alist = alist.reset_index().rename({0: 'amps'}, axis=1)
    ahlist = df.groupby('high').apply(lambda x: set(x.amp))
    ahlist = ahlist.reset_index().rename({0: 'amps'}, axis=1)
    for v1, v2 in tqdm(habitat_list):
        o1 = alist.loc[alist['general_envo_name'] == v1, 'amps'].tolist()[0]
        o2 = alist.loc[alist['general_envo_name'] == v2, 'amps'].tolist()[0]
        o = len(o1.intersection(o2))
        tests.append((v1, v2, 'observed', o))
    for v1, v2 in tqdm(hhabitat_list):
        o1 = ahlist.loc[ahlist['high'] == v1, 'amps'].tolist()[0]
        o2 = ahlist.loc[ahlist['high'] == v2, 'amps'].tolist()[0]
        o = len(o1.intersection(o2))
        testsh.append((v1, v2, 'observed', o))
    del alist, ahlist
    print('Starting permutations')
    ndf = df.groupby(['sample', 'general_envo_name']).apply(lambda x: set(x.amp))
    ndf = ndf.reset_index()
    ndf = ndf.rename({0: 'amp'}, axis=1)
    labels = ndf['general_envo_name'].tolist()
    for idx in tqdm(range(n)):
        shuffle(labels)
        ndf['general_envo_name'] = labels
        ndf['high'] = [higher_level.get(x) for x in labels]
        alist = ndf.groupby('general_envo_name').apply(lambda x: set(chain.from_iterable(x.amp)))
        alist = alist.reset_index().rename({0: 'amps'}, axis=1)
        ahlist = ndf.groupby('high').apply(lambda x: set(chain.from_iterable(x.amp)))
        ahlist = ahlist.reset_index().rename({0: 'amps'}, axis=1)
        for v1, v2 in habitat_list:
            o1 = alist.loc[alist['general_envo_name'] == v1, 'amps'].tolist()[0]
            o2 = alist.loc[alist['general_envo_name'] == v2, 'amps'].tolist()[0]
            o = len(o1.intersection(o2))
            tests.append((v1, v2, 'perm', o))
        for v1, v2 in hhabitat_list:
            o1 = ahlist.loc[ahlist['high'] == v1, 'amps'].tolist()[0]
            o2 = ahlist.loc[ahlist['high'] == v2, 'amps'].tolist()[0]
            o = len(o1.intersection(o2))
            testsh.append((v1, v2, 'perm', o))
    tests = pd.DataFrame(tests, columns=['env1', 'env2', 'test', 'overlap_AMPs'])
    testsh = pd.DataFrame(testsh, columns=['env1', 'env2', 'test', 'overlap_AMPs'])
    return (tests, testsh)


def test_perm(permutations, verbose: bool = None):
    res = []
    habitat_list = set(permutations['env1']).union(set(permutations['env2']))
    habitat_list = list(combinations(habitat_list, 2))
    for v1, v2 in habitat_list:
        xdf = permutations[(permutations.env1 == v1) & (permutations.env2 == v2)]
        if len(xdf) == 0:
            xdf = permutations[(permutations.env1 == v2) & (permutations.env2 == v1)]
        observed = xdf.loc[xdf.test == 'observed',
                           'overlap_AMPs'].tolist()[0]
        perms = xdf.loc[xdf.test == 'perm',
                        'overlap_AMPs'].tolist()
        avg, std = np.mean(perms), np.std(perms)
        _, sp = shapiro(perms)  # test normal distribution
        z_score = (observed - avg) / std
        p_value = norm.sf(abs(z_score))*2  # test result
        if verbose == None: verbose = True
        if verbose:
            print(f'''The overlap between {v1} and {v2}: {observed} c_AMPs
            Permutation test: avg = {avg}, std = {std}
            Normal distribution: Shapiro-Wilk - p = {sp}
            Overlap significance: p = {p_value}''')
        res.append((v1, v2, observed, avg, std, sp, z_score, p_value))
    res = pd.DataFrame(res,
                       columns=['habitat1',
                                'habitat2',
                                'overlap_amps',
                                'random_overlap_avg',
                                'random_overlap_std',
                                'shapiro_test',
                                'z_score',
                                'p_value'])
    return res


permutations, permutationsh = permtest(df, 1000)

#exporting results
permutations.to_csv('outputs/permutation_overlap_AMPs.tsv.gz',
                    sep='\t',
                    header=True,
                    index=None)

permutationsh.to_csv('outputs/permutation_overlap_AMPs_highlevel.tsv.gz',
                     sep='\t',
                     header=True,
                     index=None)


# In[27]:


perm_res = test_perm(permutations, False)
perm_resh = test_perm(permutationsh, False)

perm_res.to_csv('outputs/permutation_overlap_AMPs_stats.tsv.gz',
                sep='\t',
                header=True,
                index=None)

perm_resh.to_csv('outputs/permutation_overlap_AMPs_highlevel_stats.tsv.gz',
                 sep='\t',
                 header=True,
                 index=None)


# ## Testing significance of c_AMPs overlap in high level groups of habitats
#
# We tested overlap of c_AMP contents of high level groups of habitats presenting AMPs. Then, a permutation test was performed by shuffling the labels of samples and recalculating the overlap between each pair of habitats 32 times. The average and standard deviation overlap c_AMPs was calculated and the Z-score of the actual measure was computed. The p-value is then calculated using the survival function of *scipy.stats*.

# **NOTE:** Shapiro-Wilk test showed that during some permutations, the distribution was not normal. To address this problem, we suggest convert the z_scores in the final dataframe by using the [Chebyshev's inequality](), which states the probability of data distant more than N standard deviations from the average:
#
# Pr(|𝜇−𝜎| >= k*𝜎) <= 1 / k^2
#
# Although it does not reveal much about the percentile position, gives a probability that can be extended to all observations regardless the distribution.

# In[28]:


# for python it can be implemented as follows:
def Chebyshev_inequality(num_std_deviations):

    return 1 / num_std_deviations**2


# ### Creating Sample-based accumulation curves
#
# We generate sets of c_AMPs per sample, then group them by high/habitat. We then randomly select a path to continuously grow the curves. This random seed is then added consecutively to randomly selected samples and the number of AMPs is counted until end the process. The curve is computed 32 times for each environment and the average is plotted with the error for each dot.

# **Collector's curve for high level habitat groups**

# selecting data
itdf = df.groupby(['high', 'general_envo_name', 'sample'])['amp'].apply(lambda x: set(x)).reset_index()
itdf


# In[30]:


high_habitats = itdf.high.value_counts()
high_habitats = high_habitats.sort_values()
high_habitats = high_habitats.index



def collectors(df, envo: str, col: str, perms=None, step=None):
    if step == None: step = 10
    if perms == None: perms = 32
    n_itdf = df[df[col] == envo]['amp']
    f = []
    L = len(n_itdf)
    for k in tqdm(range(1, L+1, step)):
        if k != L:
            presult = 0
            for p in range(perms):
                s = n_itdf.sample(k)
                s = set(chain.from_iterable(s))
                presult += len(s)
            presult = presult / perms
        else:
            presult = len(set(chain.from_iterable(n_itdf)))
        f.append(presult)
    return (envo, f)

ccurves = []
for i in high_habitats:
    h, c = collectors(df=itdf,
                      envo=i,
                      col='high',
                      perms=32)

    ccurves.append((h, c))

avg_table = pd.DataFrame([x[1] for x in ccurves],
                         index=[x[0] for x in ccurves]).T

avg_table['samples'] = [(x*10)+1 for x in avg_table.index]
avg_table = avg_table.set_index('samples')

fig, ax = plt.subplots()
sns.lineplot(data=avg_table/1000, palette='Dark2', ax=ax)

ax.set_ylabel('c_AMPs (Thousands)')
ax.set_xlabel('Random samples')
fig.savefig('outputs/collector_curves_highenvo.svg')

# export avg_table
avg_table.to_csv('outputs/collectors_curve_highenvo.tsv.gz',
                 sep='\t',
                 header=True,
                 index=True)


# only test environments with at least 100 samples
k = itdf.general_envo_name.value_counts()
k = k.sort_values()
k = k[k >= 100].index

sitdf = itdf[itdf.general_envo_name.isin(k)]

ccurves = []
for i in k:
    h, c = collectors(df=sitdf,
                      envo=i,
                      col='general_envo_name',
                      perms=32)

    ccurves.append((h, c))


avg_table = pd.DataFrame([x[1] for x in ccurves],
                         index=[x[0] for x in ccurves]).T

avg_table['samples'] = [(x*10)+1 for x in avg_table.index]
avg_table = avg_table.set_index('samples')

fig, ax = plt.subplots()
sns.lineplot(data=avg_table/1000,
             palette='Dark2',
             ax=ax)

ax.set_ylabel('c_AMPs (Thousands)')
ax.set_xlabel('Random samples')
fig.savefig('outputs/collector_curves_generalenvo.svg')


# export avg_table
avg_table.to_csv('outputs/collectors_curve_generalenvo.tsv.gz',
                 sep='\t',
                 header=True,
                 index=True)

