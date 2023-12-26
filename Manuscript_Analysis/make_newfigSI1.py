from Bio import SeqIO
from macrel.AMP_features import fasta_features

seqs = list()
for record in SeqIO.parse('AMP.train.fasta', 'fasta'):
    name, group = record.description.split('|')
    if group == 'NAMP':
        seqs.append((name, group, str(record.seq)))

with open('NAMPs.fasta', 'w') as ofile:
    for name, group, sequence in seqs:
        ofile.write(f'>{name}\n{sequence}\n')

feat = fasta_features('NAMPs.fasta')
feat = feat.reset_index().rename({'index': 'Access'}, axis=1)
feat.to_csv('macrel_trainneg.features.tsv.gz', sep='\t', header=True, index=None)

###########################

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import combinations
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
plt.rcParams['svg.fonttype'] = 'none'

ampsphere = pd.read_table('/home/celio/Documents/updated_datasets_mar2023/data_folder/ampsphere_v2022-03.features.tsv.gz', sep='\t')
dramp = pd.read_csv('/home/celio/Documents/updated_datasets_mar2023/data_folder/dramp_v3.features.tsv.gz', sep="\t")
macrel = pd.read_csv('/home/celio/Documents/updated_datasets_mar2023/data_folder/macrel_trainpos.features.tsv.gz', sep="\t")
namps = pd.read_csv('/home/celio/Documents/updated_datasets_mar2023/data_folder/macrel_trainneg.features.tsv.gz', sep="\t")

namps['group'] = 'Macrel Negative Train'
macrel['group'] = 'Macrel Positive Train'
dramp['group'] = 'DRAMP'
ampsphere['group'] = 'AMPSphere'

df = pd.concat([ampsphere, dramp, macrel, namps])
df['length'] = df.sequence.apply(lambda x: len(x))
df.drop('sequence', axis=1, inplace=True)
df.reset_index(drop=True, inplace=True)

a = df.groupby('group').agg('mean').T
b = df.groupby('group').agg('std').T

b.to_excel('features2.xlsx', sheet_name='std', header=True, index=True)
a.to_excel('features1.xlsx', sheet_name='means', header=True, index=True)

# initializing axis array, and graph_keys
fig, axarr = plt.subplot_mosaic([['a)', 'b)', 'c)'],
                                 ['d)', 'e)', 'f)'],
                                 ['g)', 'h)', 'i)']])


graphkey = {'a)': ['length', 'Length (residues)'], 
            'b)': ['smallAA', 'Small residues'],
            'c)': ['basicAA', 'Basic residues'],
            'd)': ['pI', 'Isoelectric point'],
            'e)': ['charge', 'Charge at pH 7.0'],
            'f)': ['aindex', 'Aliphatic index'],
            'g)': ['instaindex', 'Instability index'],
            'h)': ['boman', 'Boman index'],
            'i)': ['hmoment', 'Hydrophobic moment']}

palette = {'Macrel Negative Train': '#1b9e77',
           'Macrel Positive Train': '#d95f02',
           'DRAMP': '#7570b3',
	   'AMPSphere': '#e7298a'}

for k in graphkey:
   feat, label = graphkey[k]
   for group in set(df.group):
       sns.kdeplot(ax=axarr[k],
                   data=df[df.group == group],
                   x=feat,
                   label=group,
                   color=palette[group])
   axarr[k].set_xlabel(label)
   axarr[k].set_xlim(df[feat].min(), df[feat].max())
   if k == 'a)':
       axarr[k].legend()
   if k not in ['a)', 'd)', 'g)']:
       axarr[k].set_yticks([])
       axarr[k].set_ylabel(None)

for label, ax in axarr.items():
    ax.set_title(label,
                 fontfamily='Sans Serif',
                 fontsize='large',
                 loc='left')

plt.savefig('new_figSI1.svg')
plt.close()

groups = set(df.group)
vars = list(df.columns[2:])
combs = list(combinations(groups, 2))
comparisons = []
for v in vars:
    for k, j in combs:
        u, p = mannwhitneyu(df[df.group == k][v], df[df.group == j][v])
        comparisons.append((k, j, v, u, p))
        
comparisons = pd.DataFrame(comparisons, columns=['G1', 'G2', 'Variable', 'U', 'p-value'])

_, comparisons['padj'], _, _ = multipletests(comparisons['p-value'])

comparisons.to_csv('comparison_variables.tsv', sep='\t', header=True, index=None)
