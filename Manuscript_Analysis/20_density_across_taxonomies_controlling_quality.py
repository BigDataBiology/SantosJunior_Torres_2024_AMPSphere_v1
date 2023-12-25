#!/usr/bin/env python
# coding: utf-8

# # AMPSphere v.2022-03¶
# 
# This is a notebook meant to form the set of notebooks used to analyze the data in AMPSphere and write the manuscript:
# 
# __AMPSphere: Global survey of prokaryotic antimicrobial peptides shaping microbiomes__
# 
# 
# ### Analysis of c_AMP density in different genera using only quality-controlled AMPs
# 
# Most of our data shows that species more abundandant in these environments, tend to have more AMP genes. To address this issue, we normalized the number of c_AMP genes in these genera dividing it by the assembled base pairs per genera. In this analysis, we only used AMPs passing the coordinates test to ensure the results we obtained using
# the entire dataset.

# In[1]:


# load libraries
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from Bio import Phylo
from re import sub
from tqdm import tqdm
from itertools import combinations
from matplotlib.colors import to_hex
from scipy.stats import spearmanr, mannwhitneyu
from scipy.stats import percentileofscore as posr
from statsmodels.stats.multitest import multipletests


# In[2]:


# load data
# getting amp genes per taxon
quality = pd.read_table('../data_folder/quality_assessment.tsv.xz')
quality.set_index('AMP', inplace=True)
quality = list(quality[quality.Coordinates == 'Passed'].index)

taxonomy_amps = pd.read_table('../data_folder/complete_amps_associated_taxonomy.tsv.gz')
taxonomy_amps = taxonomy_amps[taxonomy_amps.amp.isin(quality)]
taxonomy_amps = taxonomy_amps[taxonomy_amps.level.isin(['species', 'genus'])]
taxonomy_amps['fixed'] = taxonomy_amps.source.apply(lambda x: x.split()[0])
taxonomy_amps = taxonomy_amps.groupby('fixed').agg('size')
taxonomy_amps = taxonomy_amps.sort_values()

bps = pd.read_table('../data_folder/bps-per-taxon.tsv.xz', header='infer')
bps = bps[bps.level.isin(['species', 'genus'])]
bps['fixed'] = bps['name'].apply(lambda x: x.split(' ')[0])


# In[3]:


# work dfs
taxbps = bps[['fixed', 'nbps']].groupby('fixed').agg('sum')
taxbps = pd.concat([taxonomy_amps, taxbps], axis=1).fillna(0)
taxbps = taxbps.rename({0: 'amp_genes'}, axis=1)
taxbps['amps_per_Gbp'] = taxbps['amp_genes'] * 1e9 / taxbps['nbps']


# In[4]:


# plot AMP per bp per taxon
sns.scatterplot(data=taxbps, x='nbps', y='amp_genes', s=3)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Assembled base pairs per genus')
plt.ylabel('AMP genes per genus')


# In[5]:


# get statistics
a = taxbps.loc[taxbps.amp_genes > 0, 'nbps'].mean() / 1e9
b = taxbps.loc[taxbps.amp_genes > 0, 'nbps'].std() / 1e9
c = taxbps.loc[taxbps.amp_genes == 0, 'nbps'].mean() / 1e9
d = taxbps.loc[taxbps.amp_genes == 0, 'nbps'].std() / 1e9

print(f'''Statistics for the assembly bps by taxon:

- Taxa presenting AMPs
Average: {a:.2f} Gbp
Stdev.: {b:.2f} Gbp

- Taxa not presenting AMPs
Average: {c:.2f} Gbp
Stdev.: {d:.2f} Gbp''')


# In[6]:


# calculating quantiles for top taxa in AMPSphere
idf = taxbps[taxbps.amps_per_Gbp > 0]
idf = idf.sort_values(by='amps_per_Gbp')
idf = idf['amps_per_Gbp']

print('Percentile of AMP densities by genera:\n')
print('Bradyrhizobium', '\t', posr(idf, idf.loc['Bradyrhizobium'], kind='rank'))
print('Prevotella', '\t', posr(idf, idf.loc['Prevotella'], kind='rank'))
print('Pelagibacter', '\t', posr(idf, idf.loc['Pelagibacter'], kind='rank'))
print('Faecalibacterium', '\t', posr(idf, idf.loc['Faecalibacterium'], kind='rank'))
print('CAG-110', '\t', posr(idf, idf.loc['CAG-110'], kind='rank'))


# In[7]:


# add error to normalized AMPs per Gbp
p = taxbps['amp_genes'] / taxbps['nbps']
taxbps['MOE'] = np.sqrt(p * (1 - p) / taxbps['nbps'])

# our confidence level is set to 95%
# then Z=1.96
taxbps['MOE'] = 1.96 * taxbps['MOE']  
# fix proportion to Gbp
taxbps['MOE'] = taxbps['MOE'] * 1e9

taxbps['UL'] = taxbps['amps_per_Gbp'] + taxbps['MOE']
taxbps['LL'] = taxbps['amps_per_Gbp'] - taxbps['MOE']
taxbps['VAR_pct'] = taxbps['MOE'] * 100 / taxbps['amps_per_Gbp']

fdata = taxbps[taxbps.VAR_pct <= 10]


# In[8]:


sns.boxplot(x=fdata['amps_per_Gbp'],
            showfliers=False,
            color='white')

sns.swarmplot(x=fdata['amps_per_Gbp'],
              color='gray',
              s=3)

toptaxa_ampsphere = ['Prevotella', 'Bradyrhizobium', 'Pelagibacter', 'Faecalibacterium', 'CAG-110']
toptaxa_density = ['Algorimicrobium', 'TMED78', 'SFJ001', 'STGJ01', 'CAG-462']

sns.swarmplot(x=fdata.loc[toptaxa_ampsphere,
                         'amps_per_Gbp'],
              color='black',
              s=5)

sns.swarmplot(x=fdata.loc[toptaxa_density,
                         'amps_per_Gbp'],
              color='red',
              s=5)

plt.xlabel('AMP density (c_AMPs per Gbp)')


# ### Checking taxonomy groups using AMP density

# In[9]:


taxbps = taxbps.drop(['UL', 'LL'], axis=1).reset_index()
taxbps.rename({'fixed': 'genus'}, axis=1, inplace=True)


# In[10]:


# prepare data
tax = pd.read_table("../data_folder/bac120_taxonomy_r207.tsv.xz",
                    sep='\t',
                    names=['domain', 'phylum',
                           'class', 'order',
                           'family', 'genus',
                           'species'])

tax_one = tax.copy()
tax_one.drop('species', axis=1, inplace=True)
tax_one = tax_one.drop_duplicates()
tax_one['genus'] = [x[3:] for x in tax_one.genus]

tax_one = tax_one.reset_index(drop=True)
tax_one = tax_one.drop_duplicates()


# In[11]:


taxbps = taxbps.merge(on='genus', right=tax_one)

# only phyla with 100 genera or more
k = taxbps.groupby('phylum').agg('size').sort_values()
k = k[k >= 100].index

fdata = taxbps[taxbps.phylum.isin(k)]
fdata = fdata.sort_values(by='phylum')

fdata = fdata[['phylum', 'genus',
               'amp_genes', 'nbps',
               'amps_per_Gbp', 'MOE',
               'VAR_pct']]


# In[12]:


# select just genera with a maximum error of 10%
fdata = fdata[fdata.VAR_pct <= 10]
fdata


# In[13]:


# eliminate outliers using Tukey fences
q1, q3 = fdata.amps_per_Gbp.quantile([0.25, 0.75])
iqr = q3 - q1
ul, ll = q3 + (1.5*iqr), q1 - (1.5*iqr)
fdata = fdata[(fdata.amps_per_Gbp >= ll) & (fdata.amps_per_Gbp <= ul)]


# In[14]:


# determine order of phyla by their Q50
porder = fdata.groupby('phylum')['amps_per_Gbp'].quantile(0.5)
porder = porder.sort_values()
print(porder)

porder = porder.index

colors = {'p__Acidobacteriota': (174/255, 199/255, 232/255, 100/255),
          'p__Actinobacteriota': (44/255, 160/255, 44/255, 100/255),
          'p__Bacteroidota': (214/255, 39/255, 40/255, 100/255),
          'p__Firmicutes': (255/255, 187/255, 120/255, 100/255),
          'p__Firmicutes_A': (31/255, 119/255, 180/255, 100/255),
          'p__Proteobacteria': (152/255, 223/255, 138/255, 100/255),
          'p__Verrucomicrobiota': (255/255, 127/255, 14/255, 100/255),
          'other': (128/255, 128/255, 128/255, 100/255)}

fdata['class'] = [k if k in colors else 'other' for k in fdata.phylum]
fdata['class'] = fdata['class'].apply(lambda x: x.replace('p__', ''))
fdata['class'] = fdata['class'].apply(lambda x: x.replace('_', ' '))


# ## Figure S5A

# In[15]:


sns.boxplot(data=fdata,
            x='class',
            y='amps_per_Gbp',
            order=['Acidobacteriota',
                   'Actinobacteriota',
                   'Bacteroidota',
                   'Firmicutes',
                   'Firmicutes A',
                   'Proteobacteria',
                   'Verrucomicrobiota',
                   'other'],
            showfliers=False,
            color='white')

palette = {k.replace('p__', '').replace('_', ' '): v for k, v in colors.items()}

sns.swarmplot(data=fdata,
              x='class',
              y='amps_per_Gbp',
              hue='class',
              order=['Acidobacteriota',
                   'Actinobacteriota',
                   'Bacteroidota',
                   'Firmicutes',
                   'Firmicutes A',
                   'Proteobacteria',
                   'Verrucomicrobiota',
                   'other'],
              palette=palette,
              s=3)
              
plt.legend().remove()
plt.xlabel('')
plt.xticks(rotation=35)
plt.ylabel('Genus-Specific c_AMP density\n(genes per assembled Gbp)')
plt.tight_layout()
plt.savefig('phyladist_density.svg')
plt.show()


# ### Analysis of c_AMP density regarding taxonomy distribution in the tree of life

# In[16]:


# load lineages from GTDB
gtdb_tree = Phylo.read('../data_folder/bac120_r202.tre', 'newick')
taxaintree = [i.name for i in gtdb_tree.get_terminals()]

tax = pd.read_table("../data_folder/bac120_taxonomy_r207.tsv.xz",
                    sep='\t',
                    names=['domain', 'phylum',
                           'class', 'order',
                           'family', 'genus',
                           'species'])

# get common genomes from tree and lineage
taxaintree = set(tax.index).intersection(set(taxaintree))
tax = tax.loc[taxaintree, 'genus']
tax = tax.reset_index()
tax.rename({'index': 'genome'}, axis=1, inplace=True)


# In[17]:


# only phyla with 100 genera or more
k = taxbps.groupby('phylum').agg('size').sort_values()
k = k[k >= 100].index

fdata = taxbps[taxbps.phylum.isin(k)]
fdata = fdata.sort_values(by='phylum')

# only taxa with maximum 10% of error
fdata = fdata[fdata.VAR_pct <= 10]

# merge data
tax.genus = [x[3:] for x in tax.genus]
tax = tax.merge(on='genus', right=fdata)
tax

