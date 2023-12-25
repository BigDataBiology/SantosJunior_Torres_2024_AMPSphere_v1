#!/usr/bin/env python
# coding: utf-8

# # AMPSphere v.2022-03
# 
# This is a notebook meant to form the set of notebooks used to analyze the data in AMPSphere and write the manuscript:
# 
# __AMPSphere: Global survey of prokaryotic antimicrobial peptides shaping microbiomes__

# ### Analysis of abundance and diversity of c_AMPs
# 
# c_AMPs are distributed as gene variants through many species. Here, we will test:
#     
#     I. Is there multi-habitat c_AMPs?
#     II. Are multi-habitat c_AMPs clonal?

# In[1]:


import lzma
import random
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from Bio import SeqIO
from tqdm import tqdm
from scipy.stats import norm
from scipy.stats import shapiro
from scipy.stats import pearsonr, spearmanr


# In[2]:


higher_level = {'sediment' : 'other',
        'bird gut' : 'other animal',
        'cat gut' : 'non-human mammal gut',
        'insect associated' : 'other animal',
        'human urogenital tract' : 'other human',
        'dog gut' : 'non-human mammal gut',
        'fermented food' : 'anthropogenic',
        'groundwater' : 'aquatic',
        'coral associated' : 'other animal',
        'rat gut' : 'non-human mammal gut',
        'human associated' : 'other human',
        'cattle gut' : 'non-human mammal gut',
        'deer gut' : 'non-human mammal gut',
        'mouse gut' : 'non-human mammal gut',
        'river associated' : 'aquatic',
        'primate gut' : 'non-human mammal gut',
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
        'pig gut' : 'non-human mammal gut',
        'human skin' : 'other human',
        'marine' : 'aquatic',
        'soil' : 'soil/plant',
        'built environment' : 'anthropogenic',
        'human gut' : 'human gut',
        'anthropogenic': 'anthropogenic',
        'bear gut' : 'non-human mammal gut',
        'bee gut': 'other animal',
        'bat gut': 'non-human mammal gut',
        'dog associated': 'other animal',
        'cattle associated': 'other animal',
        'crustacean associated': 'other animal',
        'insect gut': 'other animal',
        'goat gut': 'non-human mammal gut', 
        'rodent gut': 'non-human mammal gut',
        'fisher gut': 'non-human mammal gut',
        'human digestive tract': 'other human',
        'coyote gut': 'non-human mammal gut',
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
        'rabbit gut': 'non-human mammal gut',
        'tunicate associated': 'other animal',
        'mussel associated': 'other animal',
        'horse gut': 'non-human mammal gut',
        'wasp gut': 'other animal',
        'guinea pig gut': 'non-human mammal gut'}


# In[3]:


# load data
data = pd.read_table('../data_folder/gmsc_amp_genes_envohr_source.tsv.gz')


# In[4]:


# eliminate genes/amps without environment (from ProGenomes)
data = data[~data.general_envo_name.isna()]

# attribut high-level habitat to genes/amps
data['high'] = data.general_envo_name.map(lambda x: higher_level.get(x, 'other'))


# In[5]:


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


# In[6]:


l = get_multihabitat(data, 'high')
l0 = get_multihabitat(data, 'low')
print(f'Multi-habitat AMPs: {len(l0)}\nMulti-high-level-habitat: {len(l)}')


# In[7]:


k = data[['amp', 'high']].drop_duplicates()['amp'].value_counts()
k = k[k>1].index

with open('multihabitat_highlevel.txt', 'wt') as out:
    for i in k: out.write(f'{i}\n')

k = data[['amp', 'general_envo_name']].drop_duplicates()['amp'].value_counts()
k = k[k>1].index

with open('multihabitat_generalenvo.txt', 'wt') as out:
    for i in k: out.write(f'{i}\n')


# ### Permutation test
# 
# We shuffle the high-level habitat annotation for the samples, and then, calculate the number of multi-habitat c_AMPs. This operation is repeated 100 times, and then we calculate the average and standard deviation of the distribution of random results. Using Shapiro-Wilk test, we check if the random distribution is normal, and if it is, we calculate the Z-score for the result obtained for AMPSphere. The Z-score is then converted into a p-value to support our conclusions.

# In[8]:


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


# In[9]:


#test high level habitats
test = []
for _ in tqdm(range(100)):
    altdf = shuffle_test(data, 'high')
    altdf = len(get_multihabitat(altdf, 'high'))
    test.append(altdf)

print(test)


# In[10]:


_, p = shapiro(test)

if p < 0.05: res = 'not-normal'
else: res = 'normal'

print(f'The Shapiro-Wilks test returned a p = {p}')
print(f'This means that the distribution is {res}')

avg, std = np.mean(test), np.std(test)

print(f'Average number of random multi-habitat AMPs - {avg}, with std = {std}')

Z = (len(l) - avg) / std

pz = norm.sf(abs(Z))

print(f'The number of multi-habitat AMPs in AMPSphere was {len(l)}')
print(f'It was {Z} * std of the random distribution')
print(f'This gives us a p-value of {pz}')


# In[17]:


#test habitats
test_low = []
for _ in tqdm(range(100)):
    altdf = shuffle_test(data, 'low')
    altdf = len(get_multihabitat(altdf, 'low'))
    test_low.append(altdf)

print(test_low)


# In[18]:


_, p = shapiro(test_low)

if p < 0.05: res = 'not-normal'
else: res = 'normal'

print(f'The Shapiro-Wilks test returned a p = {p}')
print(f'This means that the distribution is {res}')

avg, std = np.mean(test_low), np.std(test_low)

print(f'Average number of random multi-habitat AMPs - {avg}, with std = {std}')

Z = (len(l0) - avg) / std

pz = norm.sf(abs(Z))

print(f'The number of multi-habitat AMPs in AMPSphere was {len(l0)}')
print(f'It was {Z} * std of the random distribution')
print(f'This gives us a p-value of {pz}')

