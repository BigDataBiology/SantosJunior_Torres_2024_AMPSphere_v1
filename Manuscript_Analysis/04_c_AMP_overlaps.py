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
from itertools import chain, permutations


# loading data
data = pd.read_table('data/gmsc_amp_genes_envohr_source.tsv.gz',
                     sep='\t',
                     header='infer')

# filter duplicates
data = data[data.is_metagenomic == True]
data = data[['amp', 'sample', 'general_envo_name']].drop_duplicates()


# eliminating environments with less than 100 samples
df = data[['sample', 'general_envo_name']].drop_duplicates()
df = df.general_envo_name.value_counts()
df = df[df >= 100].index

df = data[data.general_envo_name.isin(df)].reset_index(drop=True)
df['high'] = df.general_envo_name.map(lambda x: higher_level.get(x))
df



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


# In[26]:


permutations, permutationsh = permtest(df, 1000)

#exporting results
permutations.to_csv('permutation_overlap_AMPs.tsv.gz',
                    sep='\t',
                    header=True,
                    index=None)

permutationsh.to_csv('permutation_overlap_AMPs_highlevel.tsv.gz',
                     sep='\t',
                     header=True,
                     index=None)


# In[27]:


perm_res = test_perm(permutations, False)
perm_resh = test_perm(permutationsh, False)

perm_res.to_csv('permutation_overlap_AMPs_stats.tsv.gz',
                sep='\t',
                header=True,
                index=None)

perm_resh.to_csv('permutation_overlap_AMPs_highlevel_stats.tsv.gz',
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

# In[29]:


# selecting data
itdf = df.groupby(['high', 'general_envo_name', 'sample'])['amp'].apply(lambda x: set(x)).reset_index()
itdf


# In[30]:


high_habitats = itdf.high.value_counts()
high_habitats = high_habitats.sort_values()
high_habitats = high_habitats.index


# In[31]:


import numpy as np
from itertools import chain


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


# In[32]:


ccurves = []
for i in high_habitats:
    h, c = collectors(df=itdf,
                      envo=i,
                      col='high',
                      perms=32)
    
    ccurves.append((h, c))


# In[33]:


avg_table = pd.DataFrame([x[1] for x in ccurves],
                         index=[x[0] for x in ccurves]).T

avg_table['samples'] = [(x*10)+1 for x in avg_table.index]
avg_table = avg_table.set_index('samples')

sns.lineplot(data=avg_table/1000, palette='Dark2')

plt.ylabel('c_AMPs (Thousands)')
plt.xlabel('Random samples')
plt.savefig('collector_curves_highenvo.svg')


# In[34]:


# export avg_table
avg_table.to_csv('collectors_curve_highenvo.tsv.gz',
                 sep='\t',
                 header=True,
                 index=True)


# In[35]:


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


# In[36]:


avg_table = pd.DataFrame([x[1] for x in ccurves],
                         index=[x[0] for x in ccurves]).T

avg_table['samples'] = [(x*10)+1 for x in avg_table.index]
avg_table = avg_table.set_index('samples')

sns.lineplot(data=avg_table/1000,
             palette='Dark2')

plt.ylabel('c_AMPs (Thousands)')
plt.xlabel('Random samples')
plt.savefig('collector_curves_generalenvo.svg')


# In[37]:


# export avg_table
avg_table.to_csv('collectors_curve_generalenvo.tsv.gz',
                 sep='\t',
                 header=True,
                 index=True)

