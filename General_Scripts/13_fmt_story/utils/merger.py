def test_groups(F, df, hypothe=None):
    from scipy.stats import mannwhitneyu
    if hypothe == None: hypothe = 'two-sided'
    groups = list(df.index)
    seen = set()
    for i in groups:
        for j in groups:
            if j not in seen:
                u, p = mannwhitneyu(df.loc[i, '#_amps'],
                                    df.loc[j, '#_amps'],
                                    alternative=hypothe)
                n1, n2 = len(df.loc[i, '#_amps']), len(df.loc[j, '#_amps'])
                if p < 0.05: print(f'{F},{i},{j},{n1},{n2},{u},{p}')
            if i not in seen: seen.add(i)


def checkdf(F, df):
    import seaborn as sns
    from matplotlib import pyplot as plt
    from collections import Counter
    test = df[['#_amps', 'outcome']]
    test = test.groupby('outcome').agg(list)
    test = test.loc[test['#_amps'].apply(lambda x: len(x)) > 10]
    outcomes = [k for k,v in Counter(df.outcome).items() if v > 10]
    sns.kdeplot(data=df[df.outcome.isin(outcomes)], x='#_amps', hue='outcome')
    plt.savefig(f'{F}_all.svg')
    plt.close()
    return test_groups(F, test)


def merge_donor(df1, df2):
    '''
    df1 - represent the outcomes table
          columns: ['species', 'fmt.id', 'fmt.type',
                    'timepoint.fmt', 'sample.donor',
                    'sample.pre', 'sample.post', 'outcome',
                    'sp_retained', 'reject_consp', 'coexist_consp.rec',
                    'coexist_consp.donor', 'engraft_consp',
                    'reject_novel', 'engraft_novel', 'sp_lost',
                    'influx_consp', 'influx_novel', 'influx_total']
          
    df2 - represent the genomes/AMPs table
          columns: ['genome', 'completeness', 'contamination', 'sample',
                    'species', 'subject_type', 'timepoint.fmt',
                    'clinical_response', 'fmt.id', 'ORFs', 'smORFs',
                    '#_amps', 'sequences']
    '''
    merge_donor = df1.copy()
    merge_donor['sample'] = merge_donor['sample.donor']
    merge_donor = merge_donor.merge(on=['species', 'sample', 'fmt.id'], right=df2) # we exclude timepoint.fmt because donors do not have this info
    merge_donor.drop('sample', axis=1, inplace=True)
    merge_donor.to_csv('merged_to_donor_samples.tsv', sep='\t', header=True, index=None)
    print('Merged to donor samples -- testing AMP effects in outcome:')
    checkdf('merge_donor', merge_donor)


def merge_recipient(df1, df2):
    '''
    df1 - represent the outcomes table
          columns: ['species', 'fmt.id', 'fmt.type',
                    'timepoint.fmt', 'sample.donor',
                    'sample.pre', 'sample.post', 'outcome',
                    'sp_retained', 'reject_consp', 'coexist_consp.rec',
                    'coexist_consp.donor', 'engraft_consp',
                    'reject_novel', 'engraft_novel', 'sp_lost',
                    'influx_consp', 'influx_novel', 'influx_total']
          
    df2 - represent the genomes/AMPs table
          columns: ['genome', 'completeness', 'contamination', 'sample',
                    'species', 'subject_type', 'timepoint.fmt',
                    'clinical_response', 'fmt.id', 'ORFs', 'smORFs',
                    '#_amps', 'sequences']
    '''
    merge_recipient = df1.copy()
    merge_recipient['sample'] = merge_recipient['sample.post']
    merge_recipient = merge_recipient.merge(on=['species', 'sample', 'fmt.id'], right=df2).drop('sample', axis=1)
    merge_recipient.to_csv('merged_to_recipient_samples_post.tsv', sep='\t', header=True, index=None)
    print('Merged to recipient samples post-FMT -- testing AMP effects in outcome:')
    checkdf('merge_recipient_post', merge_recipient)
    
    merge_recipient = df1.copy()
    merge_recipient['sample'] = merge_recipient['sample.pre']
    merge_recipient = merge_recipient.merge(on=['species', 'sample', 'fmt.id'], right=df2).drop('sample', axis=1)
    merge_recipient.to_csv('merged_to_recipient_samples_pre.tsv', sep='\t', header=True, index=None)
    print('Merged to recipient samples pre-FMT -- testing AMP effects in outcome:')
    checkdf('merge_recipient_pre', merge_recipient)
    
    
def load_data():
    import pyreadr
    import pandas as pd
    result = pyreadr.read_r('data/dat.transmission.celio.Rdata')
    result = result['dat.transmission.celio']
    amp_data = pd.read_table('AMPs_per_genomes.tsv')
    amp_data = amp_data[['genome', 'ORFs',
                         'smORFs', '#_amps', 
                         'sequences']]
    sel = pd.read_table('metadata/selected_genomes.tsv')
    amp_data = sel.merge(on='genome', right=amp_data)
    amp_data.rename({'cluster': 'species'}, axis=1, inplace=True)
    amp_data = amp_data[~amp_data['#_amps'].isna()]
    return result, amp_data


def merger():
    result, amp_data = load_data()
    merge_donor(result, amp_data)
    merge_recipient(result, amp_data)
    
