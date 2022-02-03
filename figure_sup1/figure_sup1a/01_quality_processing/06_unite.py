# unite results

import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter

print('... loading test results')
d1 = pd.read_table('antifam_list.tsv', sep='\t', header='infer')
d2 = pd.read_table('rnacode_list.tsv', sep='\t', header='infer')
d3 = pd.read_table('metaproteome_list.tsv', sep='\t', header='infer')
d4 = pd.read_table('metatranscriptome_list.tsv', sep='\t', header='infer')
d5 = pd.read_table('coordinates_list.tsv', sep='\t', header='infer')

print('... fixing some concepts')
d3 = d3.replace('yellow', 'red')
d4 = d4.replace('yellow', 'red')

print('... merging all tests')
d6 = d1.merge(on='AMP', right=d2)
d6 = d6.merge(on='AMP', right=d3)
d6 = d6.merge(on='AMP', right=d4)
d6 = d6.merge(on='AMP', right=d5)
d6.sort_values(by='AMP', inplace=True)

print('... classifying results')
d6 = d6.replace('green', 'Passed')
d6 = d6.replace('red', 'Failed')
d6 = d6.replace('yellow', 'Not tested')

print('... returning experimental evidence of peptides')
intdf = d6[['metaproteomes', 'metatranscriptomes']].replace('Passed', 1).replace('Failed', 0)
intdf = intdf.sum(axis=1)
intdf = intdf.replace(0, 'Failed')
intdf = intdf.replace(1, 'Passed').replace(2, 'Passed')

print('... counting amps per class')
expevd = pd.DataFrame.from_dict(Counter(intdf), orient='index').T
antifam = pd.DataFrame.from_dict(Counter(d6.Antifam), orient='index').T
rna = pd.DataFrame.from_dict(Counter(d6.RNAcode), orient='index').T
terminal = pd.DataFrame.from_dict(Counter(d6.coordinates), orient='index').T

print('... dataframing')
df = pd.concat([antifam, terminal, expevd, rna])
df.index = ['Antifam', 'Terminal placement', 'Experimental evidence', 'RNAcode']
df = df.fillna(0)
df = df[['Passed', 'Failed', 'Not tested']]

print('... plot')
df.plot.barh(stacked=True, cmap='Dark2')
plt.xlabel('AMP candidates')
plt.ylabel('Quality tests')
plt.tight_layout()
plt.savefig('quality.svg')

print('... save files')
d6.to_csv('quality_assessment.tsv', sep='\t', header=True, index=None)

