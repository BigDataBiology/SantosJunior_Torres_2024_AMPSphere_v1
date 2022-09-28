def test_tuple(i: tuple):
    '''
    Receives a tuple with cluster and the subset of the df
    returns the test of AMP contents as mannwhitneyU
    if at time zero there is both, receiver and donor
    '''
    from scipy.stats import mannwhitneyu as mw
    f = i[1]['#_amps'].tolist()
    if len(f) > 1:
       u, pval = mw(f[0], f[1])
       n1, n2 = len(f[0]), len(f[1])
    else:
       u, pval = ('-', '-')
       item = i[1]['subject_type'].tolist()[0]
       if item == 'donor': n1, n2 = (len(f[0]), 0)
       if item == 'recipient': n1, n2 = (0, len(f[0]))
    return (i[0], n1, n2, u, pval)


def nresponders(i: tuple):
    '''
    Receives a tuple with cluster and the subset of the df
    returns the test of AMP contents as mannwhitneyU
    '''
    from scipy.stats import mannwhitneyu as mw
    f = i[1]['#_amps'].tolist()
    if len(f) > 1:
       u, pval = mw(f[0], f[1])
       n1, n2 = len(f[0]), len(f[1])
    else:
       u, pval = ('-', '-')
       item = i[1]['clinical_response'].tolist()[0]
       if item == 'responder': n1, n2 = (len(f[0]), 0)
       if item == 'non-responder': n1, n2 = (0, len(f[0]))
    return (i[0], n1, n2, u, pval)


def corr_tuple(i: tuple):
    '''
    Receives a tuple with cluster and the subset of the df
    returns the test of correlation between timepoint and
    AMP contents as SpearmanR. If at there is no more than 1
    point, returns nan
    '''
    from scipy.stats import spearmanr
    f = i[1][['timepoint.fmt', '#_amps']]
    r, pval = spearmanr(f['timepoint.fmt'], f['#_amps'])
    cluster = i[0]
    n = i[1]['timepoint.fmt']
    n = n.drop_duplicates()
    n = len(n)
    return (i[0], n, r, pval)


def load_data():
    import pandas as pd
    data = pd.read_table('AMPs_per_genomes.tsv')
    ffilter = (~(data.ORFs.isna()) & (data.subject_type == 'donor'))
    donor_species = set(data[ffilter]['cluster'])
    ffilter = (data.cluster.isin(donor_species)) & (data.subject_type != 'donor')
    recipient_species = set(data[ffilter]['cluster'])
    df = data[~data.ORFs.isna()]
    df = df.fillna('n.a.')
    df = df[['cluster',
             'subject_type',
             'timepoint.fmt',
             'clinical_response',
             'ORFs',
             'smORFs',
             '#_amps',
             'sequences']]
    df = df[df.cluster.isin(recipient_species)]
    df = df.sort_values(by=['cluster',
                            'subject_type',
                            'timepoint.fmt',
                            'clinical_response'])
    return df
    

def testing_cluster(df):
    import pandas as pd
    from statsmodels.stats.multitest import multipletests
    test = df[df['timepoint.fmt'] == 0]
    test = test[['cluster', 'subject_type', '#_amps']]
    test = test.groupby(['cluster', 'subject_type'])
    test = test.agg(lambda x: list(x))
    test = test.reset_index()
    fdf = []
    for i in test.groupby('cluster'):
        fdf.append(test_tuple(i))
    test = pd.DataFrame(fdf,
                        columns=['cluster',
                                 'donor_genomes',
                                 'recipient_genomes',
                                 'U_statistic',
                                 'p_value'])
    test = test[test.U_statistic != '-']
    test['correct_pval'] = multipletests(test['p_value'].tolist(),
                                         alpha=0.05,
                                         method='hs',
                                         is_sorted=False,
                                         returnsorted=False)[1]
    test.to_csv('test_amps_per_cluster_t0.tsv',
                sep='\t',
                header=True,
                index=None)
    
    
def test_gain(df):
    import pandas as pd
    from statsmodels.stats.multitest import multipletests
    test = df[(df.subject_type != 'donor')]
    fdf = []
    for i in test.groupby('cluster'):
        fdf.append(corr_tuple(i))
    test = pd.DataFrame(fdf,
                        columns=['cluster',
                                 'timepoints',
                                 'Spearman_rho',
                                 'Spearman_pval'])
    test = test.dropna()
    test['correct_pval'] = multipletests(test['Spearman_pval'].tolist(),
                                         alpha=0.05,
                                         method='hs',
                                         is_sorted=False,
                                         returnsorted=False)[1]
    test.to_csv('test_amps_per_cluster_difftimes.tsv',
                sep='\t',
                header=True,
                index=None)
    

def test_responders(df):
    import pandas as pd
    from scipy.stats import mannwhitneyu as mw
    from statsmodels.stats.multitest import multipletests
    cr = {'responder': 'responder',
          'ESBL_clear': 'responder',
          'relapse': 'non-responder',
          'C_diff_clear': 'responder',
          'non-responder': 'non-responder',
          'sustained_remission': 'responder'}
    test = df[df.clinical_response != 'n.a.']
    test['clinical_response'] = [cr[x] for x in test['clinical_response']]
    test = test[['cluster', 'clinical_response', '#_amps']]
    test = test.groupby(['cluster', 'clinical_response'])['#_amps']
    test = test.agg(lambda x: list(x))	
    test = test.reset_index()
    fdf = []
    for i in test.groupby('cluster'):
        fdf.append(nresponders(i))
    test = pd.DataFrame(fdf,
                        columns=['cluster',
                                 'genomes_responder',
                                 'genomes_non_responder',
                                 'U_statistic',
                                 'p_value'])
    test = test[~(test.U_statistic == '-')]
    test['correct_pval'] = multipletests(test['p_value'].tolist(),
                                         alpha=0.05,
                                         method='hs',
                                         is_sorted=False,
                                         returnsorted=False)[1]
    test.to_csv('test_amps_per_cluster_responders.tsv',
                sep='\t',
                header=True,
                index=None) 
    
        
def mtests():
    print('... load tables')
    df = load_data()
    print('''Testing AMP contents
    if donor and recipient at T0''')
    testing_cluster(df)
    print('''Testing if AMP contents
    raise or drop at diff. times''')
    test_gain(df)
    print('''Testing AMP contents change
    between responders and non_responders''')
    test_responders(df)
    
