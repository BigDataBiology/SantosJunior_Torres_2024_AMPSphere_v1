import gzip
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO
from itertools import combinations

quality = pd.read_table('quality_assessment.tsv.xz')

for k in quality.columns:
   with open(f'passed_{k}.txt', 'w') as ofile:
       for x in quality.loc[quality[k] == 'Passed', 'AMP']:
           ofile.write(f'{x}\n')

df = pd.read_table('homologs_table.tsv')
with open('dramp_starpep.txt', 'w') as ofile:
    for x in df[(df.DRAMP + df.StarPepDB45k) > 0]['AMPSphere']:
        ofile.write(f'{x}\n')

with open('smprot.txt', 'w') as ofile:
    for x in df[df.SmProtv2 == True]['AMPSphere']:
        ofile.write(f'{x}\n')

with open('gmgc.txt', 'w') as ofile:
    for x in df[df.GMGCv1 == True]['AMPSphere']:
        ofile.write(f'{x}\n')

hq = quality.set_index('AMP')
hq = hq == 'Passed'
hq = hq.drop(['metaproteomes', 'metatranscriptomes'], axis=1)
hq = hq.sum(axis=1)
hq = hq[hq == 3].index

with open('hq_candidates.txt', 'w') as ofile:
    for x in hq:
        ofile.write(f'{x}\n')

## Uses interactivenn to create the Venn's Diagrams

copred = pd.read_table('AMP_coprediction_AMPSphere.tsv')
sns.barplot(data=((copred.set_index('SeqID') >= 0.5).sum(axis=1) - 1).value_counts().sort_index().reset_index(), x='index', y=0, color='#1b9e77')
plt.xlabel('Other AMP prediction systems than Macrel')
plt.ylabel('Number of AMPs')
plt.savefig('panel_d.svg')
plt.close()

ampsphere = dict()
for record in SeqIO.parse(gzip.open('AMPSphere_v.2022-03.faa.gz', 'rt'), 'fasta'):
    ampsphere[str(record.seq)] = record.id

gmsc = pd.read_table('GMSC10.Macrel_05.AMPs.tsv.gz')
gmsc = gmsc.Sequence.value_counts()
gmsc = gmsc.reset_index().rename({0: 'sequence'}, axis=1)
gmsc['ampsphere'] = [ampsphere.get(idx) for idx in gmsc['index']]
gmsc = gmsc.dropna()

df = df.set_index('AMPSphere')
homologs = df[df.DRAMP == True].index
df2 = gmsc[gmsc.ampsphere.isin(homologs)]

for x in df.columns:
    l = df[df[x] == True].index
    n = len(gmsc[(gmsc['ampsphere'].isin(l)) & (gmsc.Sequence <= 5)])
    rare_amps.append((x, n*100/len(l)))

nonhomologs = df[df.sum(axis=1) == 0].index
n = len(gmsc[(gmsc.ampsphere.isin(nonhomologs)) & (gmsc.Sequence <= 5)])
rare_amps.append(('Non-homologous', n*100/len(nonhomologs)))

rare_amps = pd.DataFrame(rare_amps, columns=['database', 'pct_homologs_under_5'])
rare_amps.sort_values(by='pct_homologs_under_5', inplace=True)

sns.barplot(data=rare_amps, x='database', y='pct_homologs_under_5', color='#1b9e77')
plt.xlabel("Database")
plt.ylabel("% of homologs in 5 samples or less")
plt.savefig("panel_c.svg")
plt.close()

