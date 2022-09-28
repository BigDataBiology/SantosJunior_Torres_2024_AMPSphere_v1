def load_data():
    import pandas as pd
    donor = pd.read_table('merged_to_donor_samples.tsv')
    post = pd.read_table('merged_to_recipient_samples_post.tsv')
    pre = pd.read_table('merged_to_recipient_samples_pre.tsv')
    return donor, pre, post


def test_cutoff(df, label, cutoff):
    import pandas as pd
    import numpy as np
    from scipy.stats import fisher_exact
    from collections import Counter
    lt = Counter(df[df['#_amps'] < cutoff]['outcome'])
    ge = Counter(df[df['#_amps'] >= cutoff]['outcome'])
    lt_sum = sum(lt.values())
    ge_sum = sum(ge.values())
    outcomes = set(lt.keys()).intersection(ge.keys())
    test = []
    for o in outcomes:
        a, b = lt[o], ge[o]
        c, d = lt_sum - a, ge_sum - b
        table = np.array([[a, b], [c, d]])
        oddsr, p = fisher_exact(table, alternative='two-sided')
        test.append([label, o, oddsr, p])
    return pd.DataFrame(test, columns=['group', 'outcome', 'odds_ratio', 'pvalue'])
    

def run(donor, pre, post, cutoff):  
    import pandas as pd
    from statsmodels.stats.multitest import multipletests as mtp
    d = test_cutoff(donor, 'donor', cutoff)
    pr = test_cutoff(pre, 'pre-FMT', cutoff)
    po = test_cutoff(post, 'post-FMT', cutoff)
    df = pd.concat([d, pr, po])
    _, p, _, _ = mtp(df.pvalue, is_sorted=False, returnsorted=False)
    df['pvalue'] = p
    return df[df.pvalue < 0.05]
    
    
def stderror(p, n):
    import numpy as np
    return np.sqrt(p*(1-p)/n)
    
    
def calc_uncertainty(df, cutoff):
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
    import pandas as pd
    o = df[(df.pvalue < 0.05) & (df.group == i)]['outcome']
    if len(o) > 0:
        err = calc_uncertainty(v, cutoff)
        err = err[err.outcome.isin(o)]
    else:
        err = pd.DataFrame()
    return err


def grouped_barplot(err, ofile):
    import seaborn as sns
    from matplotlib import pyplot as plt
    ax = sns.barplot(data=err, x="outcome", y="prop", hue="group")     
    x_coords = [p.get_x() + 0.5*p.get_width() for p in ax.patches]
    y_coords = [p.get_height() for p in ax.patches]
    plt.errorbar(x=x_coords, y=y_coords, yerr=err["error"], fmt="none", c= "k")
    plt.legend()
    plt.savefig(ofile)
    plt.close()


def ftests():
    import pandas as pd
    donor, pre, post = load_data()
    for cutoff in range(1, 11):
        df = run(donor, pre, post, cutoff)
        df.to_csv(f'fisher_test_{cutoff}.tsv', sep='\t', header=True, index=None)
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
    
