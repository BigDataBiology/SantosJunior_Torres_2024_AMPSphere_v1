#!/usr/bin/env python


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import environments

for k in set(environments.higher_level) ^ set(higher_level.keys()):
    print(f'{k:50}\t{environments.higher_level[k]}')


is_host_associated = {'human gut' : True,
        'soil/plant' : False,
        'aquatic' : False,
        'anthropogenic' : False,
        'other human' : True,
        'non-human mammal gut' : True,
        'other animal' : True,
        'other' : False}


plants={'lettuce', 'monocots', 'cowpea',
        'mosses', 'pitcher plant', 'maize',
        'thale cress', 'siratro', 'grapevine',
        'Norway spruce', 'black cottonwood',
        'soy', 'french bean', 'silvergrass',
        'sorghum', 'bread wheat', 'sunflower',
        'carrot', 'lodgepole pine' 'burclover',
        'cottongrass', 'switchgrass', 'eudicots',
        'agave', 'barrelclover', 'alfalfa',
        'red fir'}


def listing(x):
    return ', '.join(x.tolist())


data = pd.read_table('../data_folder/gmsc_amp_genes_envohr_source.tsv.gz')
data = data.query('is_metagenomic')
data['higher'] = data['general_envo_name'].map(environments.higher_level)
data['host_associated'] = data['higher'].map(lambda g: is_host_associated.get(g, 'NA'))


nf = data[['amp','host_associated']]
nf = nf.drop_duplicates()
nf = nf.groupby('host_associated')
nf = nf.agg('size')
print(nf)


nf = data[data.higher.isin(['soil/plant', 'aquatic'])]
nf = nf[['amp', 'higher']]
nf = nf.drop_duplicates()
nf = nf.groupby('higher')
nf = nf.agg('size')
print(nf)


nf = data.query('higher == "anthropogenic"')
nf = nf[['amp','higher']]
nf = nf.drop_duplicates()
nf = nf.groupby('higher')
nf = nf.agg('size')
print(nf)


nf = len(set(data.query('higher == "other"')['amp']))
print(f'Other environments: {nf:,}')

metadata = pd.read_table('../data_folder/metadata.tsv.gz')
metadata.rename({'sample_accession': 'sample'}, axis=1, inplace=True)
nf = metadata[['sample', 'host_common_name']]
nf = nf.merge(on='sample', right=data)
nf[['host_common_name', 'amp']].drop_duplicates().groupby('host_common_name').agg('size').sort_values()


nf['nenvo'] = [x if x not in plants else 'plant' for x in nf.host_common_name]
nf[['nenvo', 'amp']].drop_duplicates().groupby('nenvo').agg('size').sort_values()


samples = data[['sample', 'higher']]
samples = samples.drop_duplicates()
samples = samples.groupby('higher')
samples = samples.agg('size')


habitats = data[['higher', 'general_envo_name']].drop_duplicates()
habitats = habitats.groupby('higher')['general_envo_name'].apply(listing)


redamps = data.groupby('higher').agg('size')


nramps = data[['higher', 'amp']]
nramps = nramps.drop_duplicates()
nramps = nramps.groupby('higher')
nramps = nramps.agg('size')


fams = pd.read_table('../data_folder/SPHERE_v.2022-03.levels_assessment.tsv.gz')
fams
fams.rename(columns={'AMP accession': 'amp',
                    'SPHERE_fam level III': 'family'},
            inplace=True)

fams = fams[['amp', 'family']]
data = data.merge(on='amp', right=fams)
fams = fams.groupby('family').agg('size')
fams = fams[fams >= 8].index


data = data[['higher', 'family']].drop_duplicates()
famps = data.groupby('higher').agg('size')
famp_l = data[data.family.isin(fams)].groupby('higher').agg('size')


# Here, it follows the supplementary table with info about the samples, number of redundant and non-redundant AMPs, as well as the number of clusters and families each high-level habitat affiliates.


df = pd.concat([habitats,
                samples,
                redamps,
                nramps,
                famps,
                famp_l],
               axis=1)

df = df.reset_index()
df.rename(columns={'higher': 'high level environment',
                'general_envo_name': 'habitats',
                0: 'samples',
                1: 'redundant AMPs',
                2: 'non-redundant AMPs',
                3: 'AMP clusters',
                4: 'AMP families'},
          inplace=True)


# ### Information about the samples used in AMPSphere

# load data again
data = pd.read_table('../data_folder/samples-min500k-assembly-prodigal-stats.tsv.gz')
amps = pd.read_table('../data_folder/gmsc_amp_genes_envohr_source.tsv.gz')


# filter columns
namps = amps[['amp',
              'sample',
              'general_envo_name']]


# eliminate redundancy
namps = namps.drop_duplicates()
namps = namps.groupby('sample')
namps = namps.agg('size')
namps = namps.reset_index()
namps.rename(columns={0: 'amps',
                  'sample': 'sample_accession'},
                 inplace=True)


# merge splitted data
a = data.merge(on='sample_accession',
               right=namps)

b = data[~data.sample_accession.isin(namps.sample_accession)]
data = pd.concat([a, b])
data.amps = data.amps.fillna(0)

data['amps_per_assembled_Mbp'] = data.amps * 1_000_000 / data.assembly_total_length



# more data...
envo = pd.read_table('../data_folder/metadata.tsv.xz')
gen = pd.read_table('../data_folder/general_envo_names.tsv.xz')
envo = envo.merge(on=['microontology',
                      'host_scientific_name',
                      'host_tax_id'],
                  right=gen,
                  how='outer')

envo = envo[~envo.sample_accession.isna()]

envo = envo[['sample_accession',
             'geographic_location',
             'latitude',
             'longitude',
             'general_envo_name',
             'environment_material',
            ]]

envo = envo.rename({'sample_accession': 'sample'}, axis=1)

data = data.merge(on='sample_accession', right=envo)


# supp table S1
sup1 = data[['sample_accession', 'general_envo_name',
             'inserts_raw', 'assembly_total_length',
             'assembly_N50', 'prodigal_total_orfs',
             'smORFs', 'amps']].copy()

sup1.columns = ['sample', 'habitat', 'raw inserts',
                'assembled bp', 'N50', 'ORFs+smORFs',
                'smORFs', 'non-redundant AMPs']

sup1

