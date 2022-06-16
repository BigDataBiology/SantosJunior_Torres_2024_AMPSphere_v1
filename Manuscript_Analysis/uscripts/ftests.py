def load_data():
    '''
    Load merged data of Fecal Matter Transplantation results table
    '''
    import pandas as pd
    donor = pd.read_table('data/merged_to_donor_samples.tsv.gz')
    post = pd.read_table('data/merged_to_recipient_samples_post.tsv.gz')
    pre = pd.read_table('data/merged_to_recipient_samples_pre.tsv.gz')
    return donor, pre, post


def test_cutoff(df, label, cutoff):
    '''
    Calculate the Fisher's test of the proportion
    of MAGs contaning different AMP ammounts.
    
    :input:
    df - pandas dataframe
    label - test set (string)
    cutoff - number of AMPs per MAG to use as cutoff
    
    :output:
    pandas dataframe containing the columns, such as group (test set),
    outcome, odds ratio and p-value
    '''
    import pandas as pd
    import numpy as np
    from scipy.stats import fisher_exact
    from collections import Counter
    # count MAGs using an AMP cutoff
    lt = Counter(df[df['#_amps'] < cutoff]['outcome'])
    lt_sum = sum(lt.values())
    ge = Counter(df[df['#_amps'] >= cutoff]['outcome'])
    ge_sum = sum(ge.values())
    # get outcomes in common from both groups
    outcomes = set(lt.keys()).intersection(ge.keys())
    # create the 2x2 table for fisher test
    test = []
    for o in outcomes:
        a, b = lt[o], ge[o]
        c, d = lt_sum - a, ge_sum - b
        table = np.array([[a, b],
                          [c, d]])
        oddsr, p = fisher_exact(table,
                                alternative='two-sided')
        test.append([label,
                     o,
                     oddsr,
                     p])
    return pd.DataFrame(test,
                        columns=['group', 'outcome',
                                 'odds_ratio', 'pvalue'])
    

def run(donor, pre, post, cutoff):  
    '''
    Scan data using different numbers of AMP per MAGs
    and then correct the p-values using HS
    
    :input:
    Pandas data frames from merged results with 
    donor, recipients pre- and post-FMT
    
    cutoff is an integer and refers to the 
    number of AMPs per MAG to use as a cutoff for
    the Fisher's test.
    
    :output:
    returns a pandas dataframe with significant p-values
    containing the tested set, outcome, odds ratio and p-value
    '''
    import pandas as pd
    from statsmodels.stats.multitest import multipletests as mtp
    # calculate fisher test by dataset having a given cutoff
    d = test_cutoff(donor,
                    'donor',
                    cutoff)
    pr = test_cutoff(pre,
                     'pre-FMT',
                     cutoff)
    po = test_cutoff(post,
                     'post-FMT',
                     cutoff)
    df = pd.concat([d, pr, po])
    # correct p-values
    _, df['pvalue'], _, _ = mtp(df.pvalue, is_sorted=False, returnsorted=False)
    return df[df.pvalue < 0.05]
    
    
def stderror(p, n):
    import numpy as np
    return np.sqrt(p*(1-p)/n)
    
    
def calc_uncertainty(df, cutoff):
    '''
    Calculate standard error of proportions
    in different sets
    
    :input:
    df - pandas dataframe
    cutoff - number of AMPs per MAGs to divide the sets and make 
             Fisher's exact test
    
    :output:
    pandas dataframe with error related to the proportion
    and the groups it come from
    '''
    import pandas as pd
    from collections import Counter
    # select groups
    lt = Counter(df[df['#_amps'] < cutoff].outcome)
    ge = Counter(df[df['#_amps'] >= cutoff].outcome)
    sum_lt, sum_ge = sum(lt.values()), sum(ge.values())
    # calculate error and proportions g1
    lt = [(k, v/sum_lt) for k, v in lt.items()]
    lt_error = [stderror(v, sum_lt) for k, v in lt]
    lt = pd.DataFrame(lt, columns=['outcome', 'prop'])
    lt['error'] = lt_error
    lt['group'] = f'<{cutoff}'
    # calculate error and proportions g2
    ge = [(k, v/sum_ge) for k, v in ge.items()]
    ge_error = [stderror(v, sum_ge) for k, v in ge]
    ge = pd.DataFrame(ge, columns=['outcome', 'prop'])
    ge['error'] = ge_error
    ge['group'] = f'>={cutoff}'
    # sum up
    return pd.concat([lt, ge])
    

def get_info(df, i, v, cutoff):
    '''
    Clean Fisher's test result eliminating
    non-significant p-values by groups and
    outcomes
    
    :input:
    df - pandas dataframe
    i - test set (group), string
    v - test set dataframe
    cutoff - AMPs per MAGs (integer)
    '''
    import pandas as pd
    o = df[(df.pvalue < 0.05) & (df.group == i)]['outcome']
    if len(o) > 0:
        err = calc_uncertainty(v, cutoff)
        err = err[err.outcome.isin(o)]
    else:
        err = pd.DataFrame()
    return err


def grouped_barplot(err, ofile):
    '''
    Plot barplots with error by outcome and test sets
    
    :input:
    err - pandas dataframe containing error, outcome and proportion
          by groups
    :output:
    ofile - image file (svg)
    '''
    import seaborn as sns
    from matplotlib import pyplot as plt
    ax = sns.barplot(data=err, x="outcome", y="prop", hue="group")     
    x_coords = [p.get_x() + 0.5*p.get_width() for p in ax.patches]
    y_coords = [p.get_height() for p in ax.patches]
    plt.errorbar(x=x_coords, y=y_coords, yerr=err["error"], fmt="none", c="k")
    plt.legend()
    plt.savefig(ofile)
    plt.close()


def main():
    '''
    Perform Fecal Matter Transplantation testing of outcomes using AMPs per MAGs
    as cutoffs.
    '''
    import pandas as pd
    
    donor, pre, post = load_data()
    mag_amps = max([donor['#_amps'].max(),
                    pre['#_amps'].max(),
                    post['#_amps'].max()])
    mag_amps = int(mag_amps) + 1
    for cutoff in range(1, mag_amps):
        df = run(donor, pre, post, cutoff)
        df.to_csv(f'fisher_test_{cutoff}.tsv',
                  sep='\t',
                  header=True,
                  index=None)
        for i, v in zip(['donor', 'pre-FMT', 'post-FMT'],
                        [donor, pre, post]):
            err = get_info(df, i, v, cutoff)
            if len(err) > 0:
                grouped_barplot(err, f'{i}_{cutoff}_test.svg')        
    
    flist = ['fisher_test_10.tsv', 'fisher_test_2.tsv',
             'fisher_test_4.tsv', 'fisher_test_6.tsv',
             'fisher_test_8.tsv', 'fisher_test_1.tsv',
             'fisher_test_3.tsv', 'fisher_test_5.tsv',
             'fisher_test_7.tsv', 'fisher_test_9.tsv']
    
    df = pd.DataFrame()
    for i in flist:
        d = pd.read_table(i)
        amp = i.replace('.tsv', '').split('_')[-1]
        d['amp_cutoff'] = amp
        df = pd.concat([df, d])

    df.to_csv('fisher_test_summary.tsv', sep='\t', header=True, index=None)
    
    
if __name__ == '__main__':
    main()
    
