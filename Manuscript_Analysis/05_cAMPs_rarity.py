#!/usr/bin/env python
# coding: utf-8

# # AMPSphere v.2022-03
# 
# This is a notebook meant to form the set of notebooks used to analyze the data in AMPSphere and write the manuscript:
# 
# __AMPSphere: Global survey of prokaryotic antimicrobial peptides shaping microbiomes__
# 
# 
# ## Infinite genes model - Rarity
# 
# Here we will show how rare the c_AMPs from different environments was computed.

# In[1]:


import pickle
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from collections import Counter
from scipy.stats import pearsonr


# In[2]:


color_map = {'human gut' : (0.8509803921568627, 0.37254901960784315, 0.00784313725490196, 0.5),
             'soil/plant' : (0.10588235294117647, 0.6196078431372549, 0.4666666666666667, 0.5),
             'aquatic' : (0.4588235294117647, 0.4392156862745098, 0.7019607843137254, 0.5),
             'anthropogenic' : (0.9058823529411765, 0.1607843137254902, 0.5411764705882353, 0.5),
             'other human' : (0.4, 0.6509803921568628, 0.11764705882352941, 0.5),
             'mammal gut' : (0.9019607843137255, 0.6705882352941176, 0.00784313725490196, 0.5),
             'other animal' : (0.6509803921568628, 0.4627450980392157, 0.11372549019607843, 0.5),
             'other' : (0.4, 0.4, 0.4, 0.5)}


# In[3]:


def metadata():
    '''
    Import metadata
    Outputs:
    1 = dictionary with keys as samples and values as the general environments
    2 = dictionary with keys as samples and values as the high-level environments
    '''
    with open('../data_folder/envo.pkl', 'rb') as handle:
        envo = pickle.load(handle)
    return (envo['general_envo_name'], envo['high'])    


# In[4]:


def calc_rare(infile):
    print('# counting nsamples')
    df = pd.read_table(infile)
    m = df['nsamples'].mean()
    n = sum(df['nsamples'] < m)
    q50 = df['nsamples'].median()
    print(f'It was observed that the median number of samples in which a AMP is detected is {q50}')
    print(f'It was observed that the mean number of samples in which a AMP is detected is {m}')
    print(f'It was observed that {n} AMPs were detected less than or equal to the mean number of detections across all AMPSphere candidates')
    print(f'It represents {n*100/len(df):.2f}% of all AMPs')


# In[5]:


def plot_test_rare(infile, envdict, ofile):
    cutoffs = [0.01, 0.1, 1, 10] 
    envdict['total'] = sum(envdict.values())
    df = pd.read_table(infile, sep='\t', header='infer')
    df = df.set_index('AMP')
    df['total'] = df.sum(axis=1)
    headers = ['Habitat',
               'Pearson_R',
               'P-value',
               'N samples',
               'Total AMPs detected',
               '0.01%_samples',
               '0.1%_samples',
               '1%_samples',
               '10%_samples']
    with open(f'{ofile}.tsv', 'w') as handle:
        handle.write('\t'.join(headers)+'\n')
        for col in df.columns:
            if envdict[col] >= 100:
                ncs = [int(envdict[col]*c/100) for c in cutoffs]
                totamps = len(df[df[col] > 0])
                a = df.loc[df[col] > 0, col].value_counts()
                a = a.sort_index()
                a = a.reset_index()
                a.columns = ['samples', 'amps']
                a['1/k'] = 1/a['samples']
                r, p = pearsonr(a.amps, a['1/k'])
                n = []
                for c in ncs:
                    i = a.loc[a['samples'] <= c, 'amps'].sum()
                    n.append(i)
                n = [str(x) for x in n]
                vals = [col, f'{r:.3f}', f'{p:.2E}',
                        str(envdict[col]),
                        str(totamps),
                        n[0], n[1], n[2], n[3]]
                handle.write('\t'.join(vals)+'\n')
                if col in color_map:
                    sns.scatterplot(y=a['amps'],
                                    x=a['samples'],
                                    label=col,
                                    alpha=0.25,
                                    color=color_map.get(col),
                                    s=3)    
                else:
                    sns.scatterplot(y=a['amps'],
                                    x=a['samples'],
                                    label=col,
                                    alpha=0.25,
                                    s=3)    
        plt.yscale('log')
        plt.xscale('log')
        plt.xlabel('Detections\n(Number of samples)')
        plt.ylabel('c_AMPs')
        plt.tight_layout()
        plt.savefig(f'{ofile}.svg')
        plt.close()
        _ = envdict.pop('total')


# In[6]:


h, hi = metadata()
h = dict(Counter(h.values()))
hi = dict(Counter(hi.values()))    


# In[7]:


calc_rare('../data_folder/unique_map_cAMPs_nsamples.tsv.xz')


# In[8]:


plot_test_rare('../data_folder/unique_map_cAMPs_high_habs.tsv.xz',
               hi,
               'output_hi')


# In[9]:


plot_test_rare('../data_folder/unique_map_cAMPs_general_habs.tsv.xz',
               h,
               'output_gen')

