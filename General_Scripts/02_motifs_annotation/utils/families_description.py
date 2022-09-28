import pandas as pd
from collections import Counter
from itertools import chain

motifs = pd.read_table('data/AMPSphere_v.2022-03.annotation.tsv')
motifs = motifs[['id', 'motif_match']]
motifs = motifs.fillna('w_o_motif')
motifs.motif_match = motifs.motif_match.apply(lambda x: x.split('|'))

newdf = pd.DataFrame()
m = set(chain.from_iterable(motifs.motif_match))
for s in m:
    motifs[s] = [ 1 if s in x else 0 for x in motifs.motif_match ]

motifs = motifs.drop('motif_match', axis=1)
motifs.rename({'id': 'amp'}, axis=1, inplace=True)

# load families
fams = pd.read_table('data/SPHERE_v.2022-03.levels_assessment.tsv.gz')
fams = fams[['AMP accession', 'SPHERE_fam level III']]
fams = fams.rename({'AMP accession': 'amp', 'SPHERE_fam level III': 'family'}, axis=1)

motifs = motifs.merge(on='amp', right=fams)
keep = [k for k, v in Counter(motifs.family).items() if v >= 8]
motifs = motifs[motifs.family.isin(keep)]
fam_sizes = dict(Counter(motifs.family))

motifs = motifs.drop('amp', axis=1).groupby('family').agg('sum')
motifs = motifs.T
for c in motifs.columns:
    motifs[c] = motifs[c] / fam_sizes[c]

motifs = motifs.T * 100

# loading quality families
qual = pd.read_table('data/quality_families.txt', sep='\t', header='infer')
qual = qual[(qual.total >= 8)]
qual = qual[(qual.perc >= 75) | (qual.experimental_evidence == True)]
selmotifs = motifs.loc[qual.family, :]

novel_families = len(selmotifs[selmotifs.w_o_motif >= 75])

dist = []
for c in selmotifs.columns:
    n = len(selmotifs[selmotifs[c] >= 75])
    dist.append((c, n))
    
dist = pd.DataFrame(dist, columns=['motif', 'families']).sort_values(by='families')
print(dist)
