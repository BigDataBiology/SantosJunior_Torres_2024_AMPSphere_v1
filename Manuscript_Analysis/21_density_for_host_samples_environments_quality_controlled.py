#!/usr/bin/env python
# coding: utf-8

# # AMPSphere v.2022-03
# 
# This is a notebook meant to form the set of notebooks used to analyze the data in AMPSphere and write the manuscript:
# 
# __AMPSphere: Global survey of prokaryotic antimicrobial peptides shaping microbiomes__
# 

# ### Check the c_AMP density of species through host vs. non-host-associated samples and different environments controlling-quality of AMPs through the coordinates test

# In[1]:


import pandas as pd
import numpy as np
import seaborn as sns
import geopandas as gpd
from matplotlib import cm
from scipy import stats

import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none'


# In[2]:


higher_level = {'sediment' : 'other',
        'bird gut' : 'other animal',
        'cat gut' : 'mammal gut',
        'insect associated' : 'other animal',
        'human urogenital tract' : 'other human',
        'dog gut' : 'mammal gut',
        'fermented food' : 'anthropogenic',
        'groundwater' : 'aquatic',
        'coral associated' : 'other animal',
        'rat gut' : 'mammal gut',
        'human associated' : 'other human',
        'cattle gut' : 'mammal gut',
        'deer gut' : 'mammal gut',
        'mouse gut' : 'mammal gut',
        'river associated' : 'aquatic',
        'primate gut' : 'mammal gut',
        'human respiratory tract' : 'other human',
        'cattle rumen' : 'mammal gut',
        'human saliva' : 'other human',
        'activated sludge' : 'anthropogenic',
        'lake associated' : 'aquatic',
        'wastewater' : 'anthropogenic',
        'chicken gut' : 'other animal',
        'air' : 'other',
        'human mouth' : 'other human',
        'plant associated' : 'soil/plant',
        'water associated' : 'aquatic',
        'pig gut' : 'mammal gut',
        'human skin' : 'other human',
        'marine' : 'aquatic',
        'soil' : 'soil/plant',
        'built environment' : 'anthropogenic',
        'human gut' : 'human gut',
        'anthropogenic': 'anthropogenic',
        'bear gut' : 'mammal gut'}


is_host_associated = {'human gut' : True,
        'soil/plant' : False,
        'aquatic' : False,
        'anthropogenic' : False,
        'other human' : True,
        'mammal gut' : True,
        'other animal' : True,
        'other' : False}


custom_palette = {'host': (208/255, 94/255, 21/255, 100/255),
                  'non-host': (112/255, 107/255, 170/255, 100/255)}


# In[3]:


print('loading data')
meta = pd.read_table('../data_folder/metadata.tsv.xz')
syns = pd.read_table('../data_folder/general_envo_names.tsv.xz')
samples = pd.read_table('../data_folder/samples-min500k-assembly-prodigal-stats.tsv.xz', index_col=0)
gmsc = pd.read_table("../data_folder/gmsc_amp_genes_envohr_source.tsv.gz")
amp = pd.read_table('../data_folder/complete_amps_associated_taxonomy.tsv.gz')
quality = pd.read_table('../data_folder/quality_assessment.tsv.xz')


# In[4]:


print('preparing data')

quality.set_index('AMP', inplace=True)
quality = list(quality[quality.Coordinates == 'Passed'].index)

meta = meta.merge(syns[['general_envo_name',
                        'host_tax_id',
                        'microontology']],
                  on=['microontology',
                      'host_tax_id'])
                      
meta['higher'] = meta['general_envo_name'].map(lambda g: higher_level.get(g, 'other'))
meta.set_index('sample_accession', inplace=True)

amp = amp.loc[amp.amp.isin(quality), 'gmsc'].tolist()
gmsc = gmsc[(gmsc.is_metagenomic == True) & (gmsc.gmsc.isin(amp))]
gmsc = gmsc['sample'].value_counts()
samples = pd.concat([samples, gmsc], axis=1)
cols = list(samples.columns)
cols[-1] = 'amp_genes'
samples.columns = cols
del cols

for c in ['inserts_filtered', 'smORFs', 'assembly_total_length', 'amp_genes']:
    meta[c] = samples[c]

meta['smorfs_per_assembly_mbps'] = meta.eval('1_000_000 * smORFs/assembly_total_length')
meta['amp_genes_per_assembly_mbps'] = meta.eval('1_000_000 * amp_genes/assembly_total_length')

meta['is_host_associated'] = meta['general_envo_name'].map(lambda c : is_host_associated[higher_level.get(c, 'other')])
meta['is_host_associated'] = meta.is_host_associated.map(lambda i: 'host' if i else 'non-host')

assembly_total_length = meta.groupby('general_envo_name').sum()['assembly_total_length']


# In[5]:


c_general_envo_name = meta['general_envo_name'].value_counts()
sel = meta[meta.general_envo_name.map(lambda e: c_general_envo_name[e] >= 100)]
sel = sel.query('assembly_total_length > 1_000_000')
sel = sel.query('amp_genes_per_assembly_mbps > 0')
order = sel.groupby('general_envo_name').median()['amp_genes_per_assembly_mbps'].sort_values().index

sels = []
for h in order:
    cur = sel[(sel.general_envo_name == h) & (sel.amp_genes_per_assembly_mbps < 4)]
    sels.append(cur.sample(100, replace=True))

sell2000=pd.concat(sels)


# ### Figure S5B

# In[6]:


sns.boxplot(x='general_envo_name',
        y='amp_genes_per_assembly_mbps',
        order=order,
        color='white',
        showfliers=False,
        data=sel)

sns.stripplot(x='general_envo_name',
        y='amp_genes_per_assembly_mbps',
        hue='is_host_associated',
        order=order,
        data=sell2000,
        s=1.5,
        palette=custom_palette)

plt.xticks(range(32), rotation=90)
plt.xlabel('Habitat')
plt.ylabel('AMPs per assembled Mbp')
plt.tight_layout()
plt.show()


# ### Figure S5C

# In[7]:


sns.boxplot(data=sel,
           x='is_host_associated',
           y='amp_genes_per_assembly_mbps',
           order=['host', 'non-host'], 
           showfliers=False,
           color='white')

sns.swarmplot(x='is_host_associated',
              y='amp_genes_per_assembly_mbps',
              order=['host', 'non-host'],
              data=sell2000,
              s=1.5,
              palette=custom_palette)

plt.xticks(range(2), rotation=90)
plt.xlabel('Habitat type')
plt.ylabel('AMPs per assembled Mbp')
plt.tight_layout()
plt.show()

