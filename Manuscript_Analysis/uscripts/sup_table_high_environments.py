import pandas as pd

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

is_host_associated = {'human gut' : True,
        'soil/plant' : False,
        'aquatic' : False,
        'anthropogenic' : False,
        'other human' : True,
        'non-human mammal gut' : True,
        'other animal' : True,
        'other' : False}

data = pd.read_table('data/gmsc_amp_genes_envohr_source.tsv.gz')
data = data[data.is_metagenomic == True]
data['higher'] = data['general_envo_name'].map(lambda g: higher_level.get(g, 'other'))
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

nf = data[data.higher == 'anthropogenic']
nf = nf[['amp','higher']]
nf = nf.drop_duplicates()
nf = nf.groupby('higher')
nf = nf.agg('size')
print(nf)

nf = len(set(data[(data.higher == 'other')]['amp']))
print(f'Other environments: {nf}')

metadata = pd.read_table('data/metadata.tsv.gz')
metadata.rename({'sample_accession': 'sample'}, axis=1, inplace=True)
nf = metadata[['sample', 'host_common_name']]
nf = nf.merge(on='sample', right=data)
nf[['host_common_name', 'amp']].drop_duplicates().groupby('host_common_name').agg('size').sort_values()
print(nf)

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
       
nf['nenvo'] = [x if x not in plants else 'plant' for x in nf.host_common_name]
nf[['nenvo', 'amp']].drop_duplicates().groupby('nenvo').agg('size').sort_values()  

samples = data[['sample', 'higher']]
samples = samples.drop_duplicates()
samples = samples.groupby('higher')
samples = samples.agg('size')

def listing(x):
    return ', '.join(x.tolist())

habitats = data[['higher', 'general_envo_name']].drop_duplicates()
habitats = habitats.groupby('higher')['general_envo_name'].apply(lambda x: listing(x))

redamps = data.groupby('higher').agg('size')

nramps = data[['higher', 'amp']]
nramps = nramps.drop_duplicates()
nramps = nramps.groupby('higher')
nramps = nramps.agg('size')

fams = pd.read_table('data/SPHERE_v.2022-03.levels_assessment.tsv.gz')
fams = fams.rename({'AMP accession': 'amp',
                    'SPHERE_fam level III': 'family'},
                   axis=1)
                   
fams = fams[['amp', 'family']]
data = data.merge(on='amp', right=fams)
fams = fams.groupby('family').agg('size')
fams = fams[fams >= 8].index

data = data[['higher', 'family']].drop_duplicates()
famps = data.groupby('higher').agg('size')
famp_l = data[data.family.isin(fams)].groupby('higher').agg('size')

df = pd.concat([habitats,
                samples,
                redamps,
                nramps,
                famps,
                famp_l],
               axis=1)

df = df.reset_index()
df = df.rename({'higher': 'high level environment',
                'general_envo_name': 'habitats',
                0: 'samples',
                1: 'redundant AMPs',
                2: 'non-redundant AMPs',
                3: 'AMP clusters',
                4: 'AMP families'},
               axis=1)

df.to_csv('supplementary_table_S2.tsv',
          sep='\t',
          header=True,
          index=None)

