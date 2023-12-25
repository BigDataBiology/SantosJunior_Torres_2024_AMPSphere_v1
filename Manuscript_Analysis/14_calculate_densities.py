#!/usr/bin/env python
# coding: utf-8

# # AMPSphere v.2022-03
# 
# This is a notebook meant to form the set of notebooks used to analyze the data in AMPSphere and write the manuscript:
# 
# __AMPSphere: Global survey of prokaryotic antimicrobial peptides shaping microbiomes__
# 
# ### Check the c_AMP density of species through different environments
# 
# We generate sets of c_AMPs by species and sample. Then group them by environment and select those species happening in at least 2 habitats with at least 10 samples in each. The number of non-redundant c_AMPs attributed to a given species per sample in each sample is compared using Mann-Whitney U test and later the p-values are adjusted using Holm-Sidak.

# In[1]:


import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from itertools import combinations, chain
from scipy.stats import spearmanr, mannwhitneyu
from statsmodels.stats.multitest import multipletests as mt


# In[2]:


# load data
data = pd.read_table('../data_folder/complete_amps_associated_taxonomy.tsv.gz')
metadata = pd.read_table("../data_folder/metadata.tsv.xz")
metadata = metadata.rename({'sample_accession': 'sample'}, axis=1)


# In[3]:


# filter species-specific genes
data = data.loc[data.level == 'species', ['sample', 'source']]
data = data.groupby(['sample', 'source']).agg('size')
data = data.reset_index()
data = data.rename({'source': 'name',
                    0: 'amp_genes'}, axis=1)


# In[4]:


genenvo = pd.read_table('../data_folder/general_envo_names.tsv.xz')

metadata = metadata.merge(on=['microontology', 'host_scientific_name', 'host_tax_id'],
                          right=genenvo,
                          how='outer')

metadata = metadata[['sample', 'general_envo_name']]


# In[5]:


# calculate densities
sps = pd.read_table('../data_folder/bps-per-sample-per-species.tsv.xz',
                    sep='\t',
                    header='infer')

data = data.merge(on=['sample', 'name'], right=sps)

data['amp_density'] = data['amp_genes'] * 1e6 / data['nbps']

# add metadata
data = data.merge(on='sample', right=metadata)
data.general_envo_name = data.general_envo_name.replace('human saliva', 'human mouth')


# In[6]:


def tukeyfence(df):
    itdf = pd.DataFrame()
    for record in df.groupby(['name', 'general_envo_name']):
        record = record[1]
        if len(record) > 9:
            q1, q3 = record['amp_density'].quantile([0.25,0.75])
            iqr = q3-q1
            ul, ll = q3+(1.5*iqr), q1-(1.5*iqr)
            record = record[(record.amp_density <= ul) & (record.amp_density >= ll)]
            itdf = pd.concat([itdf, record])
    return itdf


# In[7]:


# remove outliers
df = tukeyfence(data)
df['p'] = df.amp_density / 1e6
df['moe'] = np.sqrt(df.p*(1-df.p)/df.nbps) * 1.96 * 1e6
df.drop('p', axis=1, inplace=True)
df['moe_pct'] = df.moe * 100 / df.amp_density
df.to_csv('species_amp_density_per_sample.tsv.gz',
          sep='\t',
          header=True,
          index=None)
df


# In[8]:


# select species present in >=10 samples per habitats
k = df.groupby(['name', 'general_envo_name']).agg('size')
k = k[k >= 10]

# select species present in more than 1 habitat
k = k.reset_index().groupby('name').apply(lambda x: x.general_envo_name.tolist())
k = k[k.apply(lambda x: len(x)) > 1]
k = k.reset_index()


# ### Supplementary Table S8

# In[9]:


# test difference in the amp_densitys among samples of
# the same species through different environments
test = []
for _, s, env in k.itertuples():
    combs = combinations(env, 2)
    for i, j in combs:
        n1 = df[(df.name == s) & (df.general_envo_name == i)]
        n2 = df[(df.name == s) & (df.general_envo_name == j)]
        u, p = mannwhitneyu(n1['amp_density'], n2['amp_density'])
        test.append((s, i, len(n1),
                     n1.amp_density.mean(), n1.amp_density.std(),
                     j, len(n2), n2.amp_density.mean(),
                     n2.amp_density.std(), p))

# format results
test = pd.DataFrame(test, columns=['species',
                                   'env1',
                                   'env1_samples',
                                   'amp_density_avg_e1',
                                   'amp_density_std_e1',
                                   'env2',
                                   'env2_samples',
                                   'amp_density_avg_e2',
                                   'amp_density_std_e2',
                                   'p_value'])

test

# adjust p-values
_, test['p_adj'], _, _ = mt(test['p_value'], method='hs')

# export
test.to_csv('species_amp_density_crossenvironment.tsv.gz', sep='\t', header=True, index=None)
test


# In[10]:


### Generating statistics

# select environments with at least 100 samples containing AMPs
k = data[['sample', 'general_envo_name']].drop_duplicates()
k = k.general_envo_name.value_counts()
k = k[k>=100].index

# select pairs containing both environments in the list
x = test[(test.env1.isin(k)) & (test.env2.isin(k))]

n1 = len(set(x.species))
print(f'It was tested {n1} species')

n2 = len(set(x.loc[test.p_adj < 0.05, 'species']))
print(f'''Only {round((n2*100/n1), 2)}% ({n2}) of them
presented significant differences in their AMP density
across different environments''')

savg = (x.amp_density_avg_e1.sum() + x.amp_density_avg_e2.sum()) / (2*len(x))

mmax, mmin = [max(x.amp_density_avg_e1.max(), x.amp_density_avg_e2.max()),
              min(x.amp_density_avg_e1.min(), x.amp_density_avg_e2.min())]

print(f'''In average, it was observed {savg:.2} c_AMP genes per assembled Mbp for each species
per sample, ranging between {mmax:.4}-{mmin:.2}.''')

