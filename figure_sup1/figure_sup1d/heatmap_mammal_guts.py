import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import cm
cm.YlOrBr

# file in > ubuntu@aws.big-data-biology.org:/share/work/Celio/files_for_figures/mundi_map/
data = pd.read_table('gmsc_amp_genes_envohr_source.tsv', sep='\t', header='infer')

# getting AMP lists for each host
human_gut = set(data[data['general envo name'].isin(['human gut'])]['amp'].tolist())
pig_gut = set(data[data['general envo name'].isin(['pig gut'])]['amp'].tolist())
chicken_gut = set(data[data['general envo name'].isin(['chicken gut'])]['amp'].tolist())
mouse_gut = set(data[data['general envo name'].isin(['mouse gut'])]['amp'].tolist())
dog_gut = set(data[data['general envo name'].isin(['dog gut'])]['amp'].tolist())
cat_gut = set(data[data['general envo name'].isin(['cat gut'])]['amp'].tolist())
bovine_gut = set(data[data['general envo name'].isin(['cattle gut'])]['amp'].tolist())

# stating hosts
posix_dic = {0: 'human_gut',
             1: 'pig_gut',
             2: 'chicken_gut',
             3: 'mouse_gut',
             4: 'cat_gut',
             5: 'dog_gut',
             6: 'bovine_gut'}


# calculating overlaps
matrix = [[],[],[],[],[],[],[]]
setlists = [human_gut, pig_gut, chicken_gut, mouse_gut, cat_gut, dog_gut, bovine_gut]

for n, i in enumerate(setlists):
    matrix[n].append(posix_dic[n])
    for m, j in enumerate(setlists):
        if n != m:
            matrix[n].append(len(i.intersection(j)))
        else:
            matrix[n].append(len(i))

# converting overlap info into a dataframe
df = pd.DataFrame(np.array(matrix)).set_index(0)
df.columns = posix_dic.values()
df = df.astype('int')
df.to_csv('mammal_gut_host_amp_overlap.tsv', sep='\t', header=True, index=True)

# converting overlap df into percent
df = df * 100 / df.max()
df = df.astype('int')
df.to_csv('mammal_gut_host_amp_perc_overlap.tsv', sep='\t', header=True, index=True)

# creat mask of zeros
mask = np.zeros_like(df)
mask[np.tril_indices_from(mask)] = True

# plot heatmap
sns.heatmap(df, annot=True, cmap='YlOrBr', mask=mask, square=True)
plt.tight_layout()
plt.savefig('heatmap_mammal_gut_host_overlap_amps2.svg')

