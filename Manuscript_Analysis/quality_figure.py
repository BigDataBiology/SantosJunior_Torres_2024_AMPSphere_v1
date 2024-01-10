import gzip
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from macrel.fasta import fasta_iter
from itertools import combinations

quality = pd.read_table('data/quality_assessment.tsv.xz', index_col=0)

quality = (quality == 'Passed')
quality.drop(
        ['metaproteomes', 'metatranscriptomes'],
        axis=1,
        inplace=True)

## Uses interactivenn to create the Venn's Diagrams

copred = pd.read_table('data/AMP_coprediction_AMPSphere.tsv.xz', index_col=0)
fig,ax = plt.subplots(figsize=(5,5))
sns.barplot(data=((copred >= 0.5).sum(axis=1) - 1).value_counts().sort_index().reset_index(),
            x='index',
            y='count',
            color='#1b9e77',
            ax=ax)

ax.set_xlabel('Other AMP prediction systems than Macrel')
ax.set_ylabel('Number of AMPs')
fig.tight_layout()
fig.savefig('figures/panel_d.svg')


ampsphere = {}
for h,seq in fasta_iter('data/AMPSphere_v.2022-03.faa.gz'):
    ampsphere[seq] = h

gmsc = pd.read_table('data/GMSC10.Macrel_05.AMPs.tsv.gz')
gmsc = gmsc.Sequence.value_counts()
gmsc = gmsc.reset_index().rename({0: 'sequence'}, axis=1)
gmsc['ampsphere'] = [ampsphere.get(idx) for idx in gmsc['index']]
gmsc = gmsc.dropna()

df = pd.read_table('homologs_table.tsv', index_col=0)
homologs = df[df.DRAMP].index
df2 = gmsc[gmsc.ampsphere.isin(homologs)]

rare_amps = []
for x in df.columns:
    l = df[df[x]].index
    n = len(gmsc[(gmsc['ampsphere'].isin(l)) & (gmsc.Sequence <= 5)])
    rare_amps.append((x, n*100/len(l)))

nonhomologs = df[df.sum(axis=1) == 0].index
n = len(gmsc[(gmsc.ampsphere.isin(nonhomologs)) & (gmsc.Sequence <= 5)])
rare_amps.append(('Non-homologous', n*100/len(nonhomologs)))

rare_amps = pd.DataFrame(rare_amps, columns=['database', 'pct_homologs_under_5'])
rare_amps.sort_values(by='pct_homologs_under_5', inplace=True)

fig,ax = plt.subplots()
sns.barplot(data=rare_amps, x='database', y='pct_homologs_under_5', color='#1b9e77', ax=ax)
ax.set_xlabel("Database")
ax.set_ylabel("% of homologs in 5 samples or less")
fig.savefig('figures/panel_c.svg')

