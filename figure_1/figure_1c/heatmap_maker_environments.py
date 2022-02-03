import pandas as pd
import numpy as np
from getsankey import genSankey
import plotly
import chart_studio.plotly as py

# loading data
# file in > ubuntu@aws.big-data-biology.org:/share/work/Celio/files_for_figures/mundi_map/
data = pd.read_table('gmsc_amp_genes_envohr_source.tsv', sep='\t', header='infer')

# creating sets per environment
non_mammal = set(data[data['general envo name'].isin(['bird gut', 'sponge associated', 'wasp gut', 'bee gut', 'planarian associated',  'crustacean gut', 'beatle gut', 'insect gut', 'mussel associated', 'termite gut', 'tunicate associated', 'bird skin', 'fish gut',  'ship worm associated', 'crustacean associated', 'annelidae associated', 'coral associated', 'mollusc associated', 'insect associated'])]['amp'].tolist())
mammal_gut = set(data[data['general envo name'].isin(['bear gut', 'coyote gut', 'pig gut', 'chicken gut', 'rodent gut', 'bat gut', 'fisher gut', 'rabbit gut', 'horse gut', 'rat gut', 'human gut', 'guinea pig gut', 'mouse gut', 'cattle gut', 'cat gut', 'deer gut', 'goat gut', 'primate gut', 'dog gut', 'dog associated'])]['amp'].tolist())
mammal_others = set(data[data['general envo name'].isin(['human respiratory tract', 'cattle associated', 'human mouth', 'human associated', 'human saliva', 'cattle rumen', 'human urogenital tract', 'goat rumen', 'human skin', 'human digestive tract'])]['amp'].tolist())
marine = set(data[data['general envo name'] == 'marine']['amp'].tolist())
soil = set(data[data['general envo name'] == 'soil']['amp'].tolist()) 
built_environment = set(data[data['general envo name'].isin(['built environment', 'anthropogenic'])]['amp'].tolist()) 
freshwater = set(data[data['general envo name'].isin(['ice associated', 'river associated', 'lake associated', 'groundwater', 'pond associated'])]['amp'].tolist())
wastewater = set(data[data['general envo name'].isin(['wastewater', 'activated sludge'])]['amp'].tolist()) 
plant_associated = set(data[data['general envo name'] == 'plant associated']['amp'].tolist()) 

# get other samples:
# data[~data['general envo name'].isin(['bird gut', 'sponge associated', 'wasp gut', 'bee gut', 'planarian associated',  'crustacean gut', 'beatle gut', 'insect gut', 'mussel associated', 'termite gut', 'tunicate associated', 'bird skin', 'fish gut',  'ship worm associated', 'crustacean associated', 'annelidae associated', 'coral associated', 'mollusc associated', 'insect associated', 'bear gut', 'coyote gut', 'pig gut', 'chicken gut', 'rodent gut', 'bat gut', 'fisher gut', 'rabbit gut', 'horse gut', 'rat gut', 'human gut', 'guinea pig gut', 'mouse gut', 'cattle gut', 'cat gut', 'deer gut', 'goat gut', 'primate gut', 'dog gut', 'dog associated', 'human respiratory tract', 'cattle associated', 'human mouth', 'human associated', 'human saliva', 'cattle rumen', 'human urogenital tract', 'goat rumen', 'human skin', 'human digestive tract', 'marine', 'soil', 'built environment', 'anthropogenic', 'ice associated', 'river associated', 'lake associated', 'groundwater', 'pond associated', 'wastewater', 'activated sludge', 'plant associated'])][['sample', 'general envo name']].dropna().drop_duplicates()['general envo name']
# non mammal samples >> data[data['general envo name'].isin(['bird gut', 'sponge associated', 'wasp gut', 'bee gut', 'planarian associated',  'crustacean gut', 'beatle gut', 'insect gut', 'mussel associated', 'termite gut', 'tunicate associated', 'bird skin', 'fish gut',  'ship worm associated', 'crustacean associated', 'annelidae associated', 'coral associated', 'mollusc associated', 'insect associated'])]['sample', 'general envo name'].dropna().drop_duplicates()['sample']
# mammal gut >> data[data['general envo name'].isin(['bear gut', 'coyote gut', 'pig gut', 'chicken gut', 'rodent gut', 'bat gut', 'fisher gut', 'rabbit gut', 'horse gut', 'rat gut', 'human gut', 'guinea pig gut', 'mouse gut', 'cattle gut', 'cat gut', 'deer gut', 'goat gut', 'primate gut', 'dog gut', 'dog associated'])][['sample', 'general envo name']].dropna().drop_duplicates()['sample']
# mammal others >> data[data['general envo name'].isin(['human respiratory tract', 'cattle associated', 'human mouth', 'human associated', 'human saliva', 'cattle rumen', 'human urogenital tract', 'goat rumen', 'human skin', 'human digestive tract'])][['sample', 'general envo name']].dropna().drop_duplicates()['sample']
# marine >> data[data['general envo name'] == 'marine'][['sample', 'general envo name']].dropna().drop_duplicates()['sample']
# soil >> data[data['general envo name'] == 'soil'][['sample', 'general envo name']].dropna().drop_duplicates()['sample']
# built_environment >> data[data['general envo name'].isin(['built environment', 'anthropogenic'])][['sample', 'general envo name']].dropna().drop_duplicates()['sample']
# freshwater >> data[data['general envo name'].isin(['ice associated', 'river associated', 'lake associated', 'groundwater', 'pond associated'])][['sample', 'general envo name']].dropna().drop_duplicates()['sample']
# wastewater >> data[data['general envo name'].isin(['wastewater', 'activated sludge'])][['sample', 'general envo name']].dropna().drop_duplicates()['sample']
# plant_associated >> data[data['general envo name'] == 'plant associated'][['sample', 'general envo name']].dropna()['sample'].drop_duplicates()

# calculate overlap
posix_dic = {0: 'non_mammal',
             1: 'mammal_gut',
             2: 'mammal_others',
             3: 'marine',
             4: 'soil',
             5: 'built_environment',
             6: 'freshwater',
             7: 'wastewater',
             8: 'plant_associated'}

matrix = [[],[],[],[],[],[],[],[],[]]
setlists = [non_mammal, mammal_gut, mammal_others, marine, soil, built_environment, freshwater, wastewater, plant_associated]
for n, i in enumerate(setlists):
    matrix[n].append(posix_dic[n])
    for m, j in enumerate(setlists):
        if n != m:
            matrix[n].append(len(i.intersection(j)))
        else:
            matrix[n].append(len(i))

# convert overlap info into a dataframe
df = pd.DataFrame(np.array(matrix)).set_index(0)
df.columns = posix_dic.values()
df = df.astype('int')
df.to_csv('env_amp_overlap.tsv', sep='\t', header=True, index=True)

# convert overlap into percent
df = df * 100 / df.max()
df.to_csv('env_amp_perc_overlap.tsv', sep='\t', header=True, index=True)

# create mask of zeros
mask = np.zeros_like(df)
mask[np.tril_indices_from(mask)] = True

# plot heatmap
sns.heatmap(df.astype('int'), annot=True, cmap="YlOrBr", mask=mask, square=True)
plt.tight_layout()
plt.savefig('heatmap_overlap_amps2.svg')

