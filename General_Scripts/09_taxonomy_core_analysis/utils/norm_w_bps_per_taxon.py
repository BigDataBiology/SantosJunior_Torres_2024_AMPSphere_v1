import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from scipy.stats import mannwhitneyu
#from scipy.stats import percentileofscore    

taxonomy_amps = pd.read_table('analysis/taxonomy_annotation.tsv', sep='\t', header='infer')
taxonomy_amps = taxonomy_amps[['fixed', 'amp']].drop_duplicates().groupby('fixed')['amp'].agg('size')

bps = pd.read_table('data/bps-per-taxon.tsv.xz', header='infer')
bps = bps[bps.level.isin(['species', 'genus'])]
bps['fixed'] = bps['name'].apply(lambda x: x.split(' ')[0])

taxbps = bps[['fixed', 'nbps']].groupby('fixed').agg('sum')
taxbps = pd.concat([taxonomy_amps, taxbps], axis=1).fillna(0)
taxbps['amps_per_Gbp'] = taxbps['amp'] * 1e9 / taxbps['nbps']

sns.scatterplot(data=taxbps, x='nbps', y='amp', s=3)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Assembled base pairs per genus')
plt.ylabel('Non-redundant AMPs per genus')
plt.savefig('analysis/figures/amps_per_Gbp_per_genus.svg')
plt.savefig('analysis/figures/amps_per_Gbp_per_genus.png', dpi=300)

taxbps = taxbps.reset_index()
taxbps.to_csv('analysis/amps_per_assembled_Gbp_per_taxon.tsv.xz', sep='\t', header=True, index=None)

taxbps['detected_amps'] = taxbps.amp.apply(lambda x: True if x > 0 else False)

taxbps.loc[taxbps.detected_amps == True, 'assembled_Gbp'].mean()
## 1.122858888780895
taxbps.loc[taxbps.detected_amps == True, 'assembled_Gbp'].std()
## 8.096021908118209
taxbps.loc[taxbps.detected_amps == False, 'assembled_Gbp'].mean()
## 0.014427576610062894
taxbps.loc[taxbps.detected_amps == False, 'assembled_Gbp'].std()
## 0.014027937842058226

spearmanr(taxbps.amp, taxbps.assembled_Gbp)
## SpearmanrResult(correlation=0.9389469735681419, pvalue=0.0)

q1, q3 = taxbps.nbps.quantile([0.25, 0.75])
lc, uc = (q1 - 1.5*(q3-q1)), (q3+1.5*(q3-q1))
ntaxb = taxbps[(taxbps.nbps < uc) & (taxbps.nbps > lc)].reset_index()

mannwhitneyu(ntaxb[ntaxb.detected_amps == True]['assembled_Gbp'], ntaxb[ntaxb.detected_amps == False]['assembled_Gbp'])
## MannwhitneyuResult(statistic=3143235.0, pvalue=5.236757517272394e-224)

fig, ax = plt.subplots()
sns.boxplot(x='detected_amps',
            y='assembled_Gbp',
            order=[True, False],
            ax=ax,
            showfliers=False,
            data=ntaxb,
            color='white',
           )
sns.swarmplot(x='detected_amps',
              y='assembled_Gbp',
              order=[True, False],
              ax=ax,
              data=ntaxb.sample(1000),
             )
ax.set_xlabel('Presence of AMPs within genus')
ax.set_ylabel('Assembed Gbp per genus')
fig.savefig('analysis/figures/Presence_absence_amps_assembled_Gbp.svg')
fig.savefig('analysis/figures/Presence_absence_amps_assembled_Gbp.png', dpi=300)

# calculating quantiles
idf = taxbps[taxbps.detected_amps == True].sort_values(by='amps_per_Gbp')['amps_per_Gbp']
sum(idf < idf.loc['Bradyrhizobium']) * 100 / len(idf)
## 55.13905683192261
sum([idf < idf.loc['Prevotella']) * 100 / len(idf)
## 21.753325272067716
sum(idf < idf.loc['Pelagibacter']) * 100 / len(idf)
## 47.04957678355502
sum(idf < idf.loc['Faecalibacterium']) * 100 / len(idf)
## 6.239419588875453
sum(idf < idf.loc['CAG-110']) * 100 / len(idf)
## 50.1813784764208

## calculating for different prevotella species
taxonomy_amps = pd.read_table('analysis/taxonomy_annotation.tsv', sep='\t', header='infer')
taxonomy_amps['source'] = taxonomy_amps['source'].str.replace('Prevotellamassilia', 'Prevotella massilia')

prev = taxonomy_amps[(taxonomy_amps.level == 'species') & (taxonomy_amps.source.str.contains('Prevotella'))]
prev['fixed'] = prev['source'].apply(lambda x: ' '.join(x.split(' ')[0:2]))
prev = prev[~prev.fixed.str.contains(' sp')]
prev = prev[['fixed', 'amp']].drop_duplicates().groupby('fixed').agg('size')

bps['name'] = bps['name'].str.replace('Prevotellamassilia', 'Prevotella massilia')
bprev = bps[(bps.level == 'species') & (bps['name'].str.contains('Prevotella'))]
bprev['fixed'] = bprev['name'].apply(lambda x: ' '.join(x.split(' ')[0:2]))
bprev = bprev[~bprev.fixed.str.contains(' sp')]
bpred = bprev.set_index('fixed').drop(['taxid', 'level', 'name'], axis=1)

bprev = bprev.reset_index()
prev = prev.reset_index()
prev = prev.merge(on='fixed', right=bprev).rename({0: 'amps'}, axis=1)
prev['amps_per_assembled_Gbp'] = prev['amps'] * 1e9 / prev['nbps']

# loading habitat info
sps = pd.read_table('data/prevotella_species_list.tsv')
sps['Species name'] = sps['Species name'].str.replace('P. ', 'Prevotella ')
sps = sps.rename({'Species name': 'fixed', 'Host(s)': 'host', 'Host sites(s)': 'habitat', 'Strain': 'strain'}, axis=1)
prev.to_csv('analysis/prevotella_species.tsv', sep='\t', header=True, index=None)

# add error to normalized AMPs per Gbp
data = pd.read_table('analysis/amps_per_assembled_Gbp_per_taxon.tsv.xz')
data['p'] = data['amp'] / data['nbps']
data['MOE'] = np.sqrt(data['p'] * data['1-p'] / data['nbps'])
data['MOE'] = 1.96 * data['MOE']
data['UL'] = (data['p'] + data['MOE'])*1e9
data['LL'] = (data['p'] - data['MOE'])*1e9
data = data.dropna()

def testcutoff(n, data):
    ndf = data[data['MOE'] / data['p'] <= n]
    a = set(ndf.sort_values(by='amps_per_Gbp').tail(10)['fixed'].tolist())
    b = set(ndf.sort_values(by='UL').tail(10)['fixed'].tolist())
    c = set(ndf.sort_values(by='LL').tail(10)['fixed'].tolist())
    p1, p2, p3 = len(a.intersection(b)) / 10, len(a.intersection(c)) / 10, len(b.intersection(c)) / 10
    print(f'''For cutoff of {n*100} % of margin of error, we obtained a conservation of:
    {p1} - between the AMP density and its UL
    {p2} - between the AMP density and its LL
    {p3} - between the UL and its LL of the AMP density
    ** the conservation of assessed in the top 10''')
    return n, len(ndf), a, b, c, p1, p2, p3

# test different cutoffs by different % of margin of error
tests = []
for n in range(1, 100):
    tests.append(testcutoff(n/100, data))

tests = pd.DataFrame(tests,
                     columns=['cutoff',
                              'dflen',
                              'top10_avg',
                              'top10_UL',
                              'top10_LL',
                              'c_avg_UL',
                              'c_avg_LL',
                              'c_UL_LL'])

tests.to_csv('screen_for_moe_pct.tsv.xz', sep='\t', header=True, index=None)

data['MOE'] = data['MOE']*1e9
data = data.drop('p', axis=1)
data.to_csv('analysis/normalized_amps_per_bp_per_taxon.tsv.xz',
            sep='\t',
            header=True,
            index=None)

