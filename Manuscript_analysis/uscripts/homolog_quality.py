import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import hypergeom

print('# load data')
dramp = pd.read_table('data/dramp_candidates.txt.gz', header=None)
smprot = pd.read_table('data/SmProt_candidates.txt.gz', header=None)
starpep = pd.read_table('data/starPepDB_candidates.txt.gz', header=None)
storfs = pd.read_table('data/STsORFs_candidates.txt.gz', header=None)
gmgc = pd.read_table('data/gmgc_candidates.txt.gz', header=None)

print('# generate 1st graph')
total_amps = 863_498
homologs = pd.concat([dramp, smprot, starpep, storfs, gmgc])
homologs = set(homologs[0])
data = [total_amps - len(homologs), len(homologs)]
classes = ['Unknown candidates', 'Homologs']
explode = [0, 0.1]

palette_color = sns.color_palette('Dark2')
plt.pie(data,
        labels=classes,
        colors=palette_color,
        explode=explode,
        autopct='%.0f%%')

plt.savefig('homologs_1.svg')

print('# generate 2nd graph')
data = [len(dramp), len(smprot), len(starpep), len(storfs), len(gmgc)]
classes = ['DRAMP v3', 'SmProt2', 'StarPepDB45k', 'STsORFs', 'GMGC v1']
data = pd.DataFrame(data, index=classes)
data = data.sort_values(by=0)
data = data.drop('STsORFs', axis=0)

q = pd.read_table('data/high_quality_candidates.txt.gz', header=None)
qq = pd.read_table('data/quality_candidates.txt.gz', header=None)
q = pd.concat([q, qq])
q = set(q[0])

qualprop = []
for k in [starpep, dramp, smprot, gmgc]:
    n1 = len(set(k[0]).intersection(q))
    qualprop.append(n1)

data['high_quality'] = qualprop
data.rename({0: 'total AMP candidates'}, axis=1, inplace=True)
data['other'] = data['total AMP candidates'] - data['high_quality']

pq = len(q) / 863498
data['prop_ampsphere'] = data['total AMP candidates'] * pq

x = [len(homologs),
     len(homologs.intersection(q)),
     len(homologs - q),
     pq*len(homologs)]
     
# data.loc['Union of databases'] = x

fig, ax = plt.subplots()
sns.barplot(ax=ax, data=data/1000, y='total AMP candidates', x=data.index, color='gray', label='Other AMP candidates')
sns.barplot(ax=ax, data=data/1000, y='high_quality', x=data.index, color='black', label='High quality AMP candidates')
ax.set_xlabel('Databases')
ax.set_ylabel('Thousands of AMP candidates')
ax.set_xticklabels(data.index, rotation=35)

for i, c in enumerate(data.prop_ampsphere):
    if i == 0:
        ax.hlines(y=c/1000, xmin=i-0.5,
                  xmax=i+0.5, linestyles='dashed',
                  color='red', label='Quality candidates in AMPSphere')
    else:
        ax.hlines(y=c/1000, xmin=i-0.5,
                  xmax=i+0.5, linestyles='dashed',
                  color='red')

plt.legend()
plt.tight_layout()
fig.savefig('bar_chart_homologs_n_qual.svg')

print('# calculate enrichment')
m = len(homologs.intersection(q))
M = len(homologs)
n = len(q)
N = 863498
p = hypergeom.sf(m, N, n, M)
e = (m*N)/(n*M)

print(f'''Enrichment of quality candidates among the annotated candidates
Proportion of quality candidates among the annotated: {m*100/M}
Proportion of quality candidates among the AMPSphere: {n*100/N}
Enrichment: {e} fold
p-value: {p}''')

