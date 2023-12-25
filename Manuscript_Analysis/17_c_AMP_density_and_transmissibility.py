#!/usr/bin/env python
# coding: utf-8

# ## Correlation between AMP density per species versus their transmissibility measured as SGB coefficients
# 
# We separated the species by their occurrences in human oral and gut habitats and then calculated their AMP densities using the total number of c_AMP genes per species and their accumulated assembled length. The density was then normalized to genes per assembled Mbp and were crossed by using the binomial nomenclature against the species previously reported by Valles-Colomer et al. (2023) that had their transmissibility measured using the SGB coefficient.

# In[1]:


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests


# We fixed the columns to be correlated against the AMP density from the Valles-Colomer et al. (2023) study.

# In[2]:


# getting the columns to correlate in the final test
cols = ['SGB_mother-offspring_transmissibility_until1y',
        'SGB_mother-offspring_transmissibility_until3y',
        'SGB_mother-offspring_transmissibility_until18y',
        'SGB_household_transmissibility',
        'SGB_intradataset_transmissibility', 
        'SGB_mother-infant_transmissibility']


# Then, we set some functions to obtain the genes related to microbial species, as well as filter them by environment and correlate them against the transmissibility data from Valles-Colomer and colleagues.

# In[3]:


def convdict(x: str):
    '''
    convdict is a function to convert a taxonomy lineage
    string into a dictionary returning the corresponding
    species or a None
    
    input:
    
    taxonomy lineage string with fields separated by | 
    with the different taxonomy levels separated by __
    
    output:
    
    string corresponding to the species or None in case
    it is not available
    ''' 
    f = dict()
    x = x.split('|')
    for j in x:
        j = j.split('__')
        f[j[0]] = j[1]
    return f.get('s', None)


def corrs(df, output:str):
    '''
    corrs - Calculates and plot the correlations between 
    species density of AMPs and their transmissibility
    
    input:
        df = dataframe of pandas containing different columns
        with transmissibility measures and one called density
        containing number of AMP genes per assembled base pairs
        also needs to contain a column named species, which is
        used as reference in correlation
        
        output = string used to save the data
    
    output:
        corrs = dataframe object with correlation measures
    '''
    corrs = []   
    for j in cols:
        if j in df.columns:
            r, p = spearmanr(df['density'],
                             df[j],
                             nan_policy='omit')
            L = len(df[j].dropna())
            corrs.append((j, r, p, L))
    corrs = pd.DataFrame(corrs,
                         columns=['Transmissibility', 'SpearmanR',
                                  'p-value', 'N'])
    for j in cols:
        if j in df.columns:
            x = df[['density', j]].dropna()
            x = x['density']
            y = df[j].dropna()
            sns.scatterplot(x=x, y=y)
            plt.xscale('log')
            plt.xlabel("AMP density per Mbp")
            plt.yscale('log')
            plt.ylabel("Transmissibility")
            plt.tight_layout()
            plt.savefig(f'{output}_{j}.svg')
            plt.close()
    df.to_csv(f'results_table_{output}.tsv',
              sep='\t',
              header=True,
              index=True)           
    corrs.to_csv(f'results_corr_{output}.tsv',
                 sep='\t',
                 header=True,
                 index=None)
    return corrs         


# First, by analyzing the human gut data:

# In[4]:


# getting species AMP density from guts
dens = pd.read_table('../data_folder/species_amp_density_per_sample.tsv.gz')
dens = dens[dens.general_envo_name.str.contains('human gut')]  # filter only human-associated envo
dens = dens.groupby('name').apply(lambda x: x[['amp_genes', 'nbps']].sum(axis=0))  # sum values for the same species
dens['density'] = 1_000_000 * dens.amp_genes / dens.nbps  # generates density per Mbp
dens = dens.reset_index(drop=False)
dens.rename({'name': 'species'}, axis=1, inplace=True)  # fix column name to merging it later

# getting species from the paper
#s4 - description of taxonomy for stool species
s4 = pd.read_excel('../data_folder/41586_2022_5620_MOESM3_ESM.xlsx',
                   sheet_name='Table S4')

s4['species'] = s4.Taxonomy.apply(lambda x: convdict(x))

s4 = s4[['SGB', 'species']]
s4 = s4[s4.species != None]

s4.species = s4.species.apply(lambda x: ' '.join(x.split('_')[0:2]))

# getting transmissibility measures for gut samples
s9 = pd.read_excel('../data_folder/41586_2022_5620_MOESM3_ESM.xlsx', sheet_name='Table S9')

s4 = s4.merge(on='SGB', right=s9)
s4 = s4.merge(on='species', right=dens)

df1 = corrs(s4, 'gut_only_gut')
df1['habitat'] = 'gut_vs_gut'


# Then, analyzing the human oral species:

# In[5]:


# getting species AMP density
dens = pd.read_table('../data_folder/species_amp_density_per_sample.tsv.gz')
dens = dens[dens.general_envo_name.str.contains('human mouth')]  # filter only human-associated envo
dens = dens.groupby('name').apply(lambda x: x[['amp_genes', 'nbps']].sum(axis=0))  # sum values for the same species
dens['density'] = 1_000_000 * dens.amp_genes / dens.nbps  # generates density per Mbp
dens = dens.reset_index(drop=False)
dens.rename({'name': 'species'}, axis=1, inplace=True)  # fix column name to merging it later

#s5 - description of taxonomy for saliva species
s5 = pd.read_excel('../data_folder/41586_2022_5620_MOESM3_ESM.xlsx',
                   sheet_name='Table S5')

s5['species'] = s5.Taxonomy.apply(lambda x: convdict(x))

s5 = s5[['SGB', 'species']]
s5 = s5[s5.species != None]

s5.species = s5.species.apply(lambda x: ' '.join(x.split('_')[0:2]))

# getting transmissibility measures for oral samples
s30 = pd.read_excel('../data_folder/41586_2022_5620_MOESM3_ESM.xlsx',
                    sheet_name='Table S30')

s5 = s5.merge(on='SGB', right=s30)
s5 = s5.merge(on='species', right=dens)

df2 = corrs(s5, 'oral_only_mouth')
df2['habitat'] = 'oral_vs_oral'


# Finally, we add the two results to a final dataframe structure and correct their p-values using Holm-Sidak method. The corrected data is then saved.

# In[6]:


df = pd.concat([df1, df2])

_, df['P_adj'], _, _ = multipletests(df['p-value'])

df.to_csv('summary_results_separated_niche.tsv',
         sep='\t',
         header=True,
         index=None)

print(df)

