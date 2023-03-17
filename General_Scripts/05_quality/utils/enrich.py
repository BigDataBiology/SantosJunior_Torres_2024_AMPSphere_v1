def filter_quality():
    '''
    Takes the quality assessment and other files to generate
    subsets of quality-controlled AMPs along with families 
    with substantial quality
    '''
    import pandas as pd

    print('load families')
    families = pd.read_table('data/SPHERE_v.2022-03.levels_assessment.tsv.gz',
                             sep='\t',
                             header='infer')
    families = families[['AMP accession','SPHERE_fam level III']]
    print('load qualities')
    quality = pd.read_table('quality_assessment.tsv', sep='\t', header='infer')
    quality.rename({'AMP': 'AMP accession'}, axis=1, inplace=True)
    print('merge')
    df = quality.merge(on='AMP accession', right=families)
    df = df.sort_values(by='SPHERE_fam level III')

    print('output hq')
    with open('high_quality_candidates.txt', 'w') as ofile:
        dhq = df[(df.metaproteomes == 'Passed')|(df.metatranscriptomes == 'Passed')]
        dhq = dhq[(dhq.Antifam == 'Passed') & (dhq.RNAcode == 'Passed') & (dhq.Coordinates == 'Passed')]
        for i in dhq['AMP accession']:
            ofile.write(i+'\n')

    print('output qc-AMPs')
    with open('quality_candidates.txt', 'w') as ofile:
        dq = df[(df.Antifam == 'Passed') & (df.RNAcode == 'Passed') & (df.Coordinates == 'Passed')]
        dq = dq[dq.metaproteomes != 'Passed']
        dq = dq[dq.metatranscriptomes != 'Passed']
        for i in dq['AMP accession']:
            ofile.write(i+'\n')

    print('output families')
    with open('quality_families.txt', 'w') as ofile:
        ofile.write('family\texperimental_evidence\tquality_candidates\ttotal\tperc\n')
        for dt in df.groupby('SPHERE_fam level III'):
            f, d = dt
            r1 = d[(d.Antifam == 'Passed') & (d.RNAcode == 'Passed') & (d.Coordinates == 'Passed')]
            r2 = d[(d.metaproteomes == 'Passed') | (d.metatranscriptomes == 'Passed')]
            r = pd.concat([r1, r2]).drop_duplicates()
            if len(r2) > 0: exp = True
            else: exp = False
            ratio = round(len(r) * 100/len(d), 2)
            ofile.write(f'{f}\t{exp}\t{len(r)}\t{len(d)}\t{ratio}\n')
    
    
def test_enrichment():
    '''
    Test the enrichment of quality candidates 
    among those AMPs with homologs in DRAMP
    '''
    import pandas as pd

    def enrichment(N, n, M, m):
        '''
        # M = Sequences of high-quality and quality in AMPsphere
        # m = Sequences of high-quality and quality in annotated_candidates hits
        # n = Sequences of annotated_candidates hits
        N = 863_498  # AMPSphere sequences
        '''
        from scipy.stats import hypergeom
        return hypergeom.sf(m, N, n, M)

    print('Load data')
    # dramp candidates are those AMPs with homologs in DRAMP
    annotated_candidates =  pd.read_table('data/DRAMP_annotation.raw.tsv.xz',
                                      header=None)

    annotated_candidates = annotated_candidates[annotated_candidates[10] <= 1e-5]
    annotated_candidates = annotated_candidates[0].drop_duplicates()

    q_candidates = pd.read_table('quality_candidates.txt',
                                 header=None)
                                 
    hq_candidates = pd.read_table('high_quality_candidates.txt', 
                                  header=None)
    
    N = 863_498  # AMPSphere sequences
    
    print('Setting enrichment analysis')
    print('Testing enrichment')
    print('Quality candidates vs. annotated_candidates/AMPSphere')
    q = pd.concat([hq_candidates, q_candidates])
    m = len(annotated_candidates[annotated_candidates.isin(q[0])])
    M = len(q)
    n = len(annotated_candidates)
    print(f'N {N}, n {n}, M {M}, m {m}')
    print(f'Ratio of qc candidates in annotated candidates: {m*100/n}')
    print(f'Ratio of qc candidates in AMPSphere: {M*100/N}')
    print(f'Enrichment: {(m*100*N)/(n*M*100)}')
    print(f'p-value: {enrichment(N, n, M, m)}')


def enrichment_analysis():
    filter_quality()
    test_enrichment()

    '''
    N 863498, n 6339, M 80213, m 1496
    Ratio of qc candidates in annotated candidates: 23.599936898564444
    Ratio of qc candidates in AMPSphere: 9.289309297763284
    Enrichment: 2.540548079742261
    p-value: 4.220811718840738e-250
    '''

