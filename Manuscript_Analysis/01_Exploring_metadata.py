#!/usr/bin/env python

# ###  Metadata exploration
#
# Exploration of metadata of metagenomes used in the AMPSphere.
#

import pandas as pd
import numpy as np
import seaborn as sns
import geopandas as gpd
import matplotlib.pyplot as plt

from matplotlib import cm
from scipy import stats
from itertools import combinations
from statsmodels.stats.multitest import multipletests
import environments

plt.rcParams['svg.fonttype'] = 'none'


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

color_map = environments.color_map


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


general = meta['general_envo_name'].value_counts()
general = general[general > 100]
general['other'] = len(meta) - general.sum()

samples = pd.read_table('../data_folder/samples-min500k-assembly-prodigal-stats.tsv.xz', index_col=0)
gmsc = pd.read_table("../data_folder/gmsc_amp_genes_envohr_source.tsv.gz")

gmsc = gmsc.query('is_metagenomic')
gmsc = gmsc.groupby('sample').agg('size')

samples = pd.concat([samples, gmsc], axis=1).rename({0: 'ampsphere_amps'}, axis=1)

meta = pd.concat([meta, samples], axis=1)

meta['smorfs_per_assembly_mbps'] = meta.eval('1_000_000 * smORFs/assembly_total_length')
meta['ampsphere_amps_per_assembly_mbps'] = meta.eval('1_000_000 * ampsphere_amps/assembly_total_length')
meta = meta[~meta.ampsphere_amps.isna()]


inserts_filtered = meta.groupby('general_envo_name').sum()['inserts_filtered']
inserts_filtered = inserts_filtered.reindex(general.index)
inserts_filtered['other'] = meta['inserts_filtered'].sum() - inserts_filtered.sum()
inserts_filtered //= 1000_000_000

smORFs = meta.groupby('general_envo_name').sum()['smORFs']
smORFs = smORFs.reindex(general.index)
smORFs['other'] = meta['smORFs'].sum() - smORFs.sum()

assembly_total_length = meta.groupby('general_envo_name').sum()['assembly_total_length']
assembly_total_length = assembly_total_length.reindex(general.index)
assembly_total_length['other'] = meta['assembly_total_length'].sum() - assembly_total_length.sum()

meta['is_host_associated'] = meta['general_envo_name'].map(lambda c : is_host_associated[higher_level.get(c, 'other')])
meta['is_host_associated'] = meta.is_host_associated.map(lambda i: 'host' if i else 'non-host')


# remove outliers
q1, q3 = meta.ampsphere_amps_per_assembly_mbps.quantile([0.25,0.75])
iqr = q3-q1
ll, ul = q1-(1.5*iqr), q3+(1.5*iqr)

sel = meta[(meta.ampsphere_amps_per_assembly_mbps <= ul) & (meta.ampsphere_amps_per_assembly_mbps >= ll)]

# exclude habitats with less than 100 samples
count_envo = sel['general_envo_name'].value_counts()
count_envo = count_envo[count_envo >= 100]

sel = sel[sel.general_envo_name.isin(count_envo.index)]  # at least 100 samples

order = sel.groupby('general_envo_name')['ampsphere_amps_per_assembly_mbps'].median().sort_values().index


# ## Figure 5A

# sampling 100 random points per habitat
sels = []
for h in order:
    cur = sel[sel.general_envo_name == h]
    sels.append(cur.sample(100))
sels = pd.concat(sels)

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

ax.set_xlabel('')
ax.set_ylabel('AMPs per assembled Mbp')


x = sel.query('higher == "anthropogenic"')
x.groupby('general_envo_name')['ampsphere_amps_per_assembly_mbps'].mean().sort_values()

u, p = stats.mannwhitneyu(sel[sel.is_host_associated == 'host']['ampsphere_amps_per_assembly_mbps'],
                          sel[sel.is_host_associated == 'non-host']['ampsphere_amps_per_assembly_mbps'])

print(f'Host vs. non-host > Mann-WhitneyU = {u}, p-value = {p}')

# ## Supplementary Figure S4

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

ax.set_xlabel('Habitat')
ax.set_ylabel('AMPs per assembled Mbp')
fig.tight_layout()


# ## Supplementary Table S1

# get supplementary info
m = meta[['ena_ers_sample_id', 'database', 'access_status', 'study', 'study_accession',
          'general_envo_name', 'higher', 'inserts_filtered', 'assembly_total_length',
          'smORFs', 'ampsphere_amps', 'is_host_associated']].copy()

m.rename(columns={'higher': 'macro_environment',
          'general_envo_name': 'micro_environment'},
         inplace=True)

m.reset_index(drop=True)

