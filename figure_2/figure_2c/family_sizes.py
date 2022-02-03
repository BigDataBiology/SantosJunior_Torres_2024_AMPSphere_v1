import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter

print('Load spheres')
# File SPHERE_v.2021-03.levels_assessment.tsv.gz is available in the AMPSphere zenodo repository
amp_fam = pd.read_table('SPHERE_v.2021-03.levels_assessment.tsv.gz', sep='\t', header='infer')
amp_fam = amp_fam[['AMP accession', 'SPHERE_fam level III']]
amp_fam = amp_fam.set_index('AMP accession')
amp_fam = amp_fam.to_dict()['SPHERE_fam level III']

print('Load genes')
# first execute codes for Fig.1a
# the file "taxonomy_annotation.tsv" will be generated and used in this script
data = pd.read_table('../figure_2a/taxonomy_annotation.tsv')
data = data[['amp', 'fixed']].drop_duplicates()


print('Convert amp into fams')
data['fam'] = [amp_fam[x] for x in df.amp]
data = data.sort_values('fam')
data = data.drop_duplicates()

print('Calculating number of genera per family')
df = pd.DataFrame.from_dict(Counter(data[['fam', 'fixed']].drop_duplicates()['fam']), orient='index')
df = df.sort_values(by=0)
df = df.rename({0: 'number of genera'}, axis=1)

print('Retrieving quality-controlled families')
# condition 1 - min of 8 seqs per family
a = Counter(pd.DataFrame.from_dict(amp_fam, orient='index')[0])
a = pd.DataFrame.from_dict(a, orient='index')
full_fam = set(a[a[0] >= 8].index)
spec_fam = set(df.index)
df.loc[full_fam.intersection(spec_fam)].hist(bins=100, grid=False)
plt.ylabel('Families with 8 members at least')
plt.xlabel('Number of different genera')
plt.savefig('full_families_numberofgenera.svg')

# condition 2 - 75% pass all tests or have experimental evidence
# File quality_families.txt available in ubuntu@aws.big-data-biology.org:/share/work/Celio/files_for_figures/quality_control
qual = pd.read_table('quality_families.txt', sep='\t', header='infer')
qual = qual[qual.total >= 8]
qual = set(qual.family)
df.loc[qual.intersection(spec_fam)].hist(bins=100, grid=False)
plt.ylabel('Quality-controlled families')
plt.xlabel('Number of different genera')
plt.savefig('qc_families_numberofgenera.svg')

df.loc[qual.intersection(spec_fam)][df >= 20].hist(bins=100, grid=False)
plt.xlabel('Number of different genera')
plt.ylabel('Quality-controlled families')
plt.savefig('qc_families_numberofgenera_ge20.svg')

print('save_file')
df = df.sort_index()
df['genera'] = [i[1].fixed.tolist() for i in data.set_index('fam').drop('amp', axis=1).groupby('fam')]
df = df.reset_index().rename({'index': 'family'}, axis=1)
df.to_csv('families_genera_counts.tsv', sep='\t', header=True, index=True)

