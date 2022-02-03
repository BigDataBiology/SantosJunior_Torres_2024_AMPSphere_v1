import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

print('Load data')
# file in > ubuntu@aws.big-data-biology.org:/share/work/Celio/files_for_figures/mundi_map/
data = pd.read_table('gmsc_amp_genes_envohr_source.tsv', sep='\t', header='infer')

# creating groups of AMPs per body site
print('Filter data')
skin = set(data[data['general envo name'] == 'human skin']['amp'].tolist())
respiratory_tract = set(data[data['general envo name'] == 'human respiratory tract']['amp'].tolist())
mouth = set(data[data['general envo name'] == 'human mouth']['amp'].tolist())
saliva = set(data[data['general envo name'] == 'human saliva']['amp'].tolist())
digestive_tract = set(data[data['general envo name'] == 'human digestive tract']['amp'].tolist())
gut = set(data[data['general envo name'] == 'human gut']['amp'].tolist())
urogenital_tract = set(data[data['general envo name'] == 'human urogenital tract']['amp'].tolist())

print('Stating variables')
posix_dic = {0: 'skin',
             1: 'respiratory_tract',
             2: 'mouth',
             3: 'saliva',
             4: 'digestive_tract',
             5: 'gut',
             6: 'urogenital_tract'}

# calculating overlaps
matrix = [[],[],[],[],[],[],[]]
setlists = [skin, respiratory_tract, mouth, saliva, digestive_tract, gut, urogenital_tract]

print('Generating matrices')
for n, i in enumerate(setlists):
    matrix[n].append(posix_dic[n])
    for m, j in enumerate(setlists):
        if n != m:
            matrix[n].append(len(i.intersection(j)))
        else:
            matrix[n].append(len(i))

# convert matrix into dataframe
print('converting matrix')
df = pd.DataFrame(np.array(matrix)).set_index(0)
df.columns = posix_dic.values()
df = df.astype('int')
df.to_csv('bodysites_amp_overlap.tsv', sep='\t', header=True, index=True)

# convert overlaps into percent
df = df * 100 / df.max()
df.to_csv('bodysites_amp_perc_overlap.tsv', sep='\t', header=True, index=True)

# creating mask of zeros
mask = np.zeros_like(df)
mask[np.tril_indices_from(mask)] = True

# plot heatmap
print('creating heatmaps')
sns.heatmap(df.astype('int'), annot=True, cmap="YlOrBr", mask=mask, square=True)
plt.tight_layout()
plt.savefig('heatmap_bodysites_amps2.svg')

# Analysis of unique AMPs per body site
# isolating the exclusive AMP sets
print('creating graphs of unique AMPs per envo')
matrix = [[],[],[],[],[],[],[]]
for n, i in enumerate(setlists):
    matrix[n].append(posix_dic[n])
    matrix[n].append(len(i))
    a = set()
    for m, j in enumerate(setlists):
        if n != m: a = a.union(j)
    matrix[n].append(len(i - a))

print('converting matrix')
df = pd.DataFrame(np.array(matrix), columns=['Body Site', 'Total', 'Unique'])
df = df.set_index('Body Site')
df = df.astype('int')
df.to_csv('bodysites_amp_counts.tsv', sep='\t', header=True, index=True)

# plotting exclusive genes sets
np.log10(df).plot.bar()
plt.ylabel('Log10(# AMPs)')
plt.tight_layout()
plt.savefig('amps_per_bodysite.svg')

