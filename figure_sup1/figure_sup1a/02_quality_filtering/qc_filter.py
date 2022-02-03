import pandas as pd
from collections import Counter

# File SPHERE_v.2021-03.levels_assessment.tsv.gz available in the AMPSphere Zenodo repository
families = pd.read_table('SPHERE_v.2021-03.levels_assessment.tsv.gz', sep='\t', header='infer')
families = families[['AMP accession','SPHERE_fam level III']]

# File quality_assessment.tsv produced by the 01_quality_processing set of scripts
quality = pd.read_table('quality_assessment.tsv', sep='\t', header='infer')
quality.rename({'AMP': 'AMP accession'}, axis=1, inplace=True)

df = quality.merge(on='AMP accession', right=families)
df = df.sort_values(by='SPHERE_fam level III')

with open('high_quality_candidates.txt', 'w') as ofile:
    dhq = df[(df.metaproteomes == 'Passed')|(df.metatranscriptomes == 'Passed')]
    dhq = dhq[(dhq.Antifam == 'Passed') & (dhq.RNAcode == 'Passed') & (dhq.coordinates == 'Passed')]
    for i in dhq['AMP accession']:
        ofile.write(i+'\n')


with open('quality_candidates.txt', 'w') as ofile:
    dq = df[(df.Antifam == 'Passed') & (df.RNAcode == 'Passed') & (df.coordinates == 'Passed')]
    dq = dq[dq.metaproteomes != 'Passed']
    dq = dq[dq.metatranscriptomes != 'Passed']
    for i in dq['AMP accession']:
        ofile.write(i+'\n')

with open('quality_families.txt', 'w') as ofile:
    ofile.write('family\texperimental_evidence\tquality_candidates\ttotal\tperc\n')
    for dt in df.groupby('SPHERE_fam level III'):
        f, d = dt
        if len(d) > 1:
            r1 = d[(d.Antifam == 'Passed') & (d.RNAcode == 'Passed') & (d.coordinates == 'Passed')]
            r2 = d[(d.metaproteomes == 'Passed') | (d.metatranscriptomes == 'Passed')]
            r = pd.concat([r1, r2]).drop_duplicates()
            if len(r2) > 0: exp = True
            else: exp = False
        else:
            r = d[(d.metaproteomes == 'Passed') | (d.metatranscriptomes == 'Passed')]
            if len(r) > 0: exp = True
            else: exp = False            
        ratio = round(len(r) * 100/len(d), 2)
        if ratio > 75.0: ofile.write(f'{f}\t{exp}\t{len(r)}\t{len(d)}\t{ratio}\n')
    
