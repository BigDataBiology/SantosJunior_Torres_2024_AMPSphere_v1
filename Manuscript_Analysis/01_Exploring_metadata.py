#!/usr/bin/env python
# coding: utf-8

# # AMPSphere manuscript
# 
# This is a notebook meant to form the set of notebooks used to analyze the data in AMPSphere and write the manuscript:
# 
# __AMPSphere: Global survey of prokaryotic antimicrobial peptides shaping microbiomes__
# 
# Figures generated in this script formed the panel shown in Figure 1A, 5A, S4 and Table S1.
# 
# 
# ###  Metadata exploration
# 
# Exploration of metadata of metagenomes used in the AMPSphere.
# 

# In[1]:


# load libraries
import pandas as pd
import numpy as np
import seaborn as sns
import geopandas as gpd
import matplotlib.pyplot as plt

from matplotlib import cm
from scipy import stats
from itertools import combinations
from statsmodels.stats.multitest import multipletests

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
        'cattle rumen' : 'other animal',
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
        'bear gut' : 'mammal gut',
        'pond associated': 'aquatic',
        'bee gut': 'other animal',
        'bat gut': 'mammal gut',
        'dog associated': 'other animal',
        'cattle associated': 'other animal',
        'crustacean associated': 'other animal',
        'insect gut': 'other animal',
        'goat gut': 'mammal gut', 
        'rodent gut': 'mammal gut',
        'fisher gut': 'mammal gut',
        'human digestive tract': 'other human',
        'coyote gut': 'mammal gut',
        'planarian associated': 'other animal',
        'sponge associated': 'other animal',
        'goat rumen': 'other animal',
        'crustacean gut': 'other animal',
        'annelidae associated': 'other animal',
        'bird skin': 'other animal',
        'beatle gut': 'other animal',
        'termite gut': 'other animal', 
        'fish gut': 'other animal',
        'mollusc associated': 'other animal',
        'ship worm associated': 'other animal',
        'rabbit gut': 'mammal gut',
        'tunicate associated': 'other animal',
        'mussel associated': 'other animal',
        'horse gut': 'mammal gut',
        'wasp gut': 'other animal',
        'guinea pig gut': 'mammal gut'}
        
is_host_associated = {'human gut' : True,
        'soil/plant' : False,
        'aquatic' : False,
        'anthropogenic' : False,
        'other human' : True,
        'mammal gut' : True,
        'other animal' : True,
        'other' : False}

color_map = {'human gut' : (0.8509803921568627, 0.37254901960784315, 0.00784313725490196),
        'soil/plant' : (0.10588235294117647, 0.6196078431372549, 0.4666666666666667),
        'aquatic' : (0.4588235294117647, 0.4392156862745098, 0.7019607843137254),
        'anthropogenic' : (0.9058823529411765, 0.1607843137254902, 0.5411764705882353),
        'other human' : (0.4, 0.6509803921568628, 0.11764705882352941),
        'mammal gut' : (0.9019607843137255, 0.6705882352941176, 0.00784313725490196),
        'other animal' : (0.6509803921568628, 0.4627450980392157, 0.11372549019607843),
        'other' : (0.4, 0.4, 0.4)}


# In[3]:


# loading data
meta = pd.read_table('../data_folder/metadata.tsv.xz')
syns = pd.read_table('../data_folder/general_envo_names.tsv.xz')

# formatting data
meta = meta.merge(syns[['general_envo_name',
                        'host_tax_id',
                        'microontology']],
                  on=['microontology', 'host_tax_id'])

meta['higher'] = meta['general_envo_name'].map(lambda g: higher_level.get(g, 'other'))
meta.set_index('sample_accession', inplace=True)
nmeta = meta[['higher', 'latitude', 'longitude']].drop_duplicates()


# ## Figure 1A

# In[4]:


fig, ax = plt.subplots(figsize=(8, 6))
countries = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
countries.plot(color="gainsboro", ax=ax)

for hab,c in color_map.items():
    sel = nmeta.query('higher == @hab')
    sel.plot(x="longitude", y="latitude", kind="scatter",
        c=[c for _ in range(len(sel))], label=hab,
        colormap="YlOrRd", s=4.5, ax=ax)

sns.despine(fig, trim=True)
ax.set_xlabel(None)
ax.set_ylabel(None)
ax.set_xticks([])
ax.set_yticks([])
ax.tick_params(left=False,
               labelleft=False,
               bottom=False,
               labelbottom=False)
#fig.show()
fig.tight_layout()
ax.legend(loc=1)


# In[5]:


general = meta['general_envo_name'].value_counts()
general = general[general > 100]
general['other'] = len(meta) - general.sum()


# In[6]:


samples = pd.read_table('../data_folder/samples-min500k-assembly-prodigal-stats.tsv.xz', index_col=0)
gmsc = pd.read_table("../data_folder/gmsc_amp_genes_envohr_source.tsv.gz")

gmsc = gmsc[gmsc.is_metagenomic == True]
gmsc = gmsc.groupby('sample').agg('size')

samples = pd.concat([samples, gmsc], axis=1).rename({0: 'ampsphere_amps'}, axis=1)

meta = pd.concat([meta, samples], axis=1)

meta['smorfs_per_assembly_mbps'] = meta.eval('1_000_000 * smORFs/assembly_total_length')
meta['ampsphere_amps_per_assembly_mbps'] = meta.eval('1_000_000 * ampsphere_amps/assembly_total_length')
meta = meta[~meta.ampsphere_amps.isna()]


# In[7]:


inserts_filtered = meta.groupby('general_envo_name').sum()['inserts_filtered']
inserts_filtered = inserts_filtered.reindex(general.index)
inserts_filtered['other'] = meta['inserts_filtered'].sum() - inserts_filtered.sum()
inserts_filtered //= 1000_000_000

smORFs = meta.groupby('general_envo_name').sum()['smORFs']
smORFs = smORFs.reindex(general.index)
smORFs['other'] = meta['smORFs'].sum() - smORFs.sum()


# In[8]:


assembly_total_length = meta.groupby('general_envo_name').sum()['assembly_total_length']
assembly_total_length = assembly_total_length.reindex(general.index)
assembly_total_length['other'] = meta['assembly_total_length'].sum() - assembly_total_length.sum()

meta['is_host_associated'] = meta['general_envo_name'].map(lambda c : is_host_associated[higher_level.get(c, 'other')])
meta['is_host_associated'] = meta.is_host_associated.map(lambda i: 'host' if i else 'non-host')


# In[9]:


# remove outliers
q1, q3 = meta.ampsphere_amps_per_assembly_mbps.quantile([0.25,0.75])
iqr = q3-q1
ll, ul = q1-(1.5*iqr), q3+(1.5*iqr)

sel = meta[(meta.ampsphere_amps_per_assembly_mbps <= ul) & (meta.ampsphere_amps_per_assembly_mbps >= ll)]


# In[10]:


# exclude habitats with less than 100 samples
count_envo = sel['general_envo_name'].value_counts()
count_envo = count_envo[count_envo >= 100]

sel = sel[sel.general_envo_name.isin(count_envo.index)]  # at least 100 samples

order = sel.groupby('general_envo_name').median()['ampsphere_amps_per_assembly_mbps'].sort_values().index


# ## Figure 5A

# In[11]:


# sampling 100 random points per habitat
sels = []
for h in order:
    cur = sel[sel.general_envo_name == h]
    sels.append(cur.sample(100))

sels=pd.concat(sels)


# In[12]:


fig,ax = plt.subplots()
ax.set_xlim(-1,2)

sns.boxplot(x='is_host_associated',
        y='ampsphere_amps_per_assembly_mbps',
        order=['host', 'non-host'],
        ax=ax,
        showfliers=False,
        data=sel.sample(2000),
        color='white',
        )

sns.swarmplot(x='is_host_associated',
        y='ampsphere_amps_per_assembly_mbps',
        order=['host', 'non-host'],
        ax=ax,
        data=sel.sample(2000),
        s=3,
        )

plt.xlabel('')
plt.ylabel('AMPs per assembled Mbp')


# In[13]:


x = sel.query('higher == "anthropogenic"')
x.groupby('general_envo_name').mean()['ampsphere_amps_per_assembly_mbps'].sort_values()

u, p = stats.mannwhitneyu(sel[sel.is_host_associated == 'host']['ampsphere_amps_per_assembly_mbps'],
                          sel[sel.is_host_associated == 'non-host']['ampsphere_amps_per_assembly_mbps'])

print(f'Host vs. non-host > Mann-WhitneyU = {u}, p-value = {p}')                          


# ## Supplementary Figure S4

# In[14]:


fig,ax = plt.subplots()

sns.boxplot(x='general_envo_name',
        y='ampsphere_amps_per_assembly_mbps',
        order=order,
        ax=ax,
        color='white',
        showfliers=False,
        data=sels)

sns.swarmplot(x='general_envo_name',
        y='ampsphere_amps_per_assembly_mbps',
        hue='is_host_associated',
        order=order,
        ax=ax,
        data=sels,
        s=1.0)

for x in ax.get_xticklabels():
    x.set_rotation(90)

plt.xlabel('Habitat')
plt.ylabel('AMPs per assembled Mbp')
plt.tight_layout()


# ## Supplementary Table S1

# In[15]:


# get supplementary info
m = meta[['ena_ers_sample_id', 'database', 'access_status', 'study', 'study_accession',
          'general_envo_name', 'higher', 'inserts_filtered', 'assembly_total_length',
          'smORFs', 'ampsphere_amps', 'is_host_associated']].copy()

m.rename({'higher': 'macro_environment',
          'general_envo_name': 'micro_environment'},
         axis=1,
         inplace=True)
          
m.reset_index(drop=True)

