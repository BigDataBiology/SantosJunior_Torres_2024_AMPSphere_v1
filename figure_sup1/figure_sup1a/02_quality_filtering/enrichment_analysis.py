import pandas as pd

print('Load data')
# File dramp_candidates.txt available in the results of Figure 1b
# annotated_candidates = pd.read_table('dramp_candidates.txt', sep='\t', header=None)
# File dramp_candidates.txt available in the results of Figure 1b
annotated_candidates = pd.read_table('../../../figure_1/figure_1b/annotated_candidates.txt', sep='\t', header=None)

# Files quality_candidates.txt and high_quality_candidates.txt are available in ubuntu@aws.big-data-biology.org:/share/work/Celio/files_for_figures/quality_control/
q_candidates = pd.read_table('quality_candidates.txt', sep='\t', header=None)
hq_candidates = pd.read_table('high_quality_candidates.txt', sep='\t', header=None)

print('Setting enrichment analysis')
# M = Sequences of high-quality and quality in AMPsphere
# m = Sequences of high-quality and quality in annotated_candidates hits
# n = Sequences of annotated_candidates hits
N = 863_498  # AMPSphere sequences

def enrichment(N, n, M, m):
    from scipy.stats import hypergeom
    return hypergeom.sf(m, N, n, M)
        
print('Testing enrichment')
print('1. high-quality candidates vs. annotated_candidates/AMPSphere')
m = len(annotated_candidates.merge(on=0, right=hq_candidates))
M = len(hq_candidates)
n = len(annotated_candidates)
print(f'N {N}, n {n}, M {M}, m {m}')
print(f'Ratio of hq candidates in annotated candidates: {m*100/n}')
print(f'Ratio of hq candidates in AMPSphere: {M*100/N}')
print(f'Enrichment: {(m*100*N)/(n*M*100)}')
print(f'p-value: {enrichment(N, n, M, m)}')

print('2. quality candidates vs. annotated_candidates/AMPSphere')
q = pd.concat([hq_candidates, q_candidates])
m = len(annotated_candidates.merge(on=0, right=q))
M = len(q)
n = len(annotated_candidates)
print(f'N {N}, n {n}, M {M}, m {m}')
print(f'Ratio of qc candidates in annotated candidates: {m*100/n}')
print(f'Ratio of qc candidates in AMPSphere: {M*100/N}')
print(f'Enrichment: {(m*100*N)/(n*M*100)}')
print(f'p-value: {enrichment(N, n, M, m)}')

'''
Load data
Setting enrichment analysis
Testing enrichment
1. high-quality candidates vs. annotated_candidates/AMPSphere
N 863498, n 73774, M 8961, m 3121
Ratio of hq candidates in annotated candidates: 4.230487705695774
Ratio of hq candidates in AMPSphere: 1.037755733076394
Enrichment: 4.076573678041389
p-value: 0.0
2. quality candidates vs. annotated_candidates/AMPSphere
N 863498, n 73774, M 30196, m 10427
Ratio of qc candidates in annotated candidates: 14.133705641553934
Ratio of qc candidates in AMPSphere: 3.49693919383716
Enrichment: 4.041736175013425
p-value: 0.0
'''

