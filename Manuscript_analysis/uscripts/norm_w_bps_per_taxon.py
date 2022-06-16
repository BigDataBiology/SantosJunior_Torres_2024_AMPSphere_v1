import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from scipy.stats import mannwhitneyu
from scipy.stats import percentileofscore as posr
from nbps_taxonomy import taxanalysis

print('# load data')
# getting amp genes per taxon
taxonomy_amps = pd.read_table('data/complete_amps_associated_taxonomy.tsv.gz')
taxonomy_amps = taxonomy_amps[taxonomy_amps.level.isin(['species', 'genus'])]
taxonomy_amps['source'] = taxonomy_amps['source'].str.replace('Prevotellamassilia', 'Prevotella massilia')
taxonomy_amps['fixed'] = taxonomy_amps.source.apply(lambda x: x.split()[0])
taxonomy_amps = taxonomy_amps.groupby('fixed').agg('size')
taxonomy_amps = taxonomy_amps.sort_values()

bps = pd.read_table('data/bps-per-taxon.tsv.xz', header='infer')
bps = bps[bps.level.isin(['species', 'genus'])]
bps['name'] = bps['name'].str.replace('Prevotellamassilia', 'Prevotella massilia')
bps['fixed'] = bps['name'].apply(lambda x: x.split(' ')[0])

print('# work dfs')
taxbps = bps[['fixed', 'nbps']].groupby('fixed').agg('sum')
taxbps = pd.concat([taxonomy_amps, taxbps], axis=1).fillna(0)
taxbps = taxbps.rename({0: 'amp_genes'}, axis=1)
taxbps['amps_per_Gbp'] = taxbps['amp_genes'] * 1e9 / taxbps['nbps']

print('# plot AMP per bp per taxon')
sns.scatterplot(data=taxbps, x='nbps', y='amp_genes', s=3)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Assembled base pairs per genus')
plt.ylabel('AMP genes per genus')
plt.savefig('amps_per_Gbp_per_genus.svg')
plt.savefig('amps_per_Gbp_per_genus.png', dpi=300)

print('# export data')
taxbps = taxbps.reset_index()
taxbps.to_csv('amps_per_assembled_Gbp_per_taxon.tsv.xz',
              sep='\t',
              header=True,
              index=None)

print('# get statistics')
a = taxbps.loc[taxbps.amp_genes > 0, 'nbps'].mean() / 1e9
b = taxbps.loc[taxbps.amp_genes > 0, 'nbps'].std() / 1e9
c = taxbps.loc[taxbps.amp_genes == 0, 'nbps'].mean() / 1e9
d = taxbps.loc[taxbps.amp_genes == 0, 'nbps'].std() / 1e9

print(f'''Statistics for the assembly bps by taxon:
- Taxa presenting AMPs
Average: {a} Gbp
Stdev.: {b} Gbp
- Taxa not presenting AMPs
Average: {c}
Stdev.: {d}''')
del a, b, c, d
 
print('# calculating quantiles for top taxa in AMPSphere')
idf = taxbps[taxbps.amps_per_Gbp > 0]
idf = idf.sort_values(by='amps_per_Gbp')
idf = idf.set_index('fixed')['amps_per_Gbp']
print('Bradyrhizobium', posr(idf, idf.loc['Bradyrhizobium'], kind='rank'))
print('Prevotella', posr(idf, idf.loc['Prevotella'], kind='rank'))
print('Pelagibacter', posr(idf, idf.loc['Pelagibacter'], kind='rank'))
print('Faecalibacterium', posr(idf, idf.loc['Faecalibacterium'], kind='rank'))
print('CAG-110', posr(idf, idf.loc['CAG-110'], kind='rank'))

print('# add error to normalized AMPs per Gbp')
taxbps['p'] = taxbps['amp_genes'] / taxbps['nbps']
taxbps['MOE'] = np.sqrt(taxbps['p'] * (1 - taxbps['p']) / taxbps['nbps'])
taxbps['MOE'] = 1.96 * taxbps['MOE'] * 1e9
taxbps['UL'] = taxbps['amps_per_Gbp'] + taxbps['MOE']
taxbps['LL'] = taxbps['amps_per_Gbp'] - taxbps['MOE']
taxbps = taxbps.drop(['p'], axis=1)
taxbps['VAR_pct'] = taxbps['MOE'] * 100 / taxbps['amps_per_Gbp']
taxbps.to_csv('normalized_amps_per_bp_per_taxon.tsv.xz',
            sep='\t',
            header=True,
            index=None)

def testcutoff(n, data):
    '''
    n - represents the max. percent error (0 to 100)
    '''
    ndf = data[data.VAR_pct <= n]
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
    tests.append(testcutoff(n, taxbps))

print('# uniting last results')
tests = pd.DataFrame(tests,
                     columns=['cutoff',
                              'dflen',
                              'top10_avg',
                              'top10_UL',
                              'top10_LL',
                              'c_avg_UL',
                              'c_avg_LL',
                              'c_UL_LL'])


print('# exporting')
tests.to_csv('screen_for_moe_pct.tsv.xz', sep='\t', header=True, index=None)

## calculating for different prevotella species
taxonomy_amps = pd.read_table('data/taxonomy_annotation.tsv.gz', sep='\t', header='infer')
taxonomy_amps['source'] = taxonomy_amps['source'].str.replace('Prevotellamassilia', 'Prevotella massilia')

prev = taxonomy_amps[(taxonomy_amps.level == 'species') & (taxonomy_amps.source.str.contains('Prevotella'))]
prev.loc[:, 'fixed'] = prev['source'].apply(lambda x: ' '.join(x.split(' ')[0:2]))
prev = prev[~prev.fixed.str.contains(' sp')]
prev = prev[['fixed', 'amp']].drop_duplicates().groupby('fixed').agg('size')

bprev = bps[(bps.level == 'species') & (bps['name'].str.contains('Prevotella'))]
bprev.loc[:, 'fixed'] = bprev['name'].apply(lambda x: ' '.join(x.split(' ')[0:2]))
bprev = bprev[~bprev.fixed.str.contains(' sp')]
bpred = bprev.set_index('fixed').drop(['taxid', 'level', 'name'], axis=1)

bprev = bprev.reset_index()
prev = prev.reset_index()
prev = prev.merge(on='fixed', right=bprev).rename({0: 'amps'}, axis=1)
prev['amps_per_assembled_Gbp'] = prev['amps'] * 1e9 / prev['nbps']

# loading habitat info
sps = pd.read_table('data/prevotella_species_list.tsv.gz')
sps = sps.rename({'Species name': 'fixed',
                  'Host(s)': 'host',
                  'Host sites(s)': 'habitat',
                  'Strain': 'strain'},
                 axis=1,
                 )
                 
sps.fixed = sps.fixed.str.replace('P. ', 'Prevotella ', regex=False)
sps = sps.merge(on='fixed', right=prev)[['name',
                                         'host',
                                         'Host site(s)',
                                         'strain',
                                         'amps',
                                         'nbps',
                                         'amps_per_assembled_Gbp']]
print(sps)
sps.to_csv('prevotella_species.tsv.gz', sep='\t', header=True, index=None)

print('# analyzing by taxonomy')
taxanalysis()

