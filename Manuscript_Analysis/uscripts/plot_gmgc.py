def main():
    print('# loading data')
    print('-- it can take a while')
    s, fpart = load_data()
    print('# plot GMGC homologs')
#    gmgc_plot(s, fpart)
    gmgc_plot2(s, fpart)


def load_data():
    import pandas as pd
    from collections import Counter
    from math import log
   
    # loading data
    gmgc = pd.DataFrame()
    for record in pd.read_table('data/result_gmgc.m8.xz',
                                header=None,
                                chunksize=1_000_000):
        gmgc = gmgc.append(record)

    gmgc.columns = ['query', 'target', 'evalue',
                    'gapopen', 'pident', 'nident', 
                    'qstart', 'qend', 'qlen',
                    'tstart', 'tend', 'tlen',
                    'alnlen', 'raw', 'bits',
                    'cigar', 'qseq', 'tseq',
                    'qheader', 'theader', 'qaln',
                    'taln', 'qframe', 'tframe',
                    'mismatch', 'qcov', 'tcov']

    print('# creating some measures')
    gmgc['pct_start'] = gmgc['tstart'] * 100 / gmgc['tlen']
    gmgc['Log(E-value)'] = [log(x, 10) for x in gmgc.evalue]
    gmgc['func'] = [x.split('.')[-1] for x in gmgc.target]

    print('# getting functions')
    f = Counter(gmgc[['query',
                      'func']].drop_duplicates()['func'])
    f = pd.DataFrame.from_dict(f,
                               orient='index').sort_values(by=0)
    fpart = f[f[0] > 400]
    fpart.loc['Other'] = f[f[0] < 400].sum()[0]
    fpart = fpart / fpart.sum()
    fpart = fpart[0].apply(lambda x: log(x, 10))
    fpart = fpart.reset_index().rename({'index': 'functions',
                                        0: 'Log(counts)'},
                                       axis=1)
    fpart = fpart.replace('UNKNOWN',
                          'Unknown')
    
    print('# separating some data')
    tmp = gmgc[['Log(E-value)', 'pident',
                'pct_start', 'qlen',
                'tlen', 'qcov',
                'tcov']]

    print('# sampling')
    s = tmp.sample(10_000).copy()
    s = s.reset_index(drop=True)
    s['Identity'] = [classify(x) for x in s.pident]

    print('# histogram of pct_start')
    import seaborn as sns
    import matplotlib.pyplot as plt
    df = pd.DataFrame()
    for r in gmgc.groupby('query'):
        r = r[1].sort_values()
        r = r.head(1)
        df = pd.concat([df, r])
    
    sns.histplot(data=df, x='pct_start', bins=100)
    plt.xlabel('Match start (% of target length)')
    plt.ylabel('Counts')
    plt.savefig('pct_start.svg')
    
    f = len(df[(df.pct_start <= 25)|(df.pct_start >= 75)])*100/len(df)
    print(f'{f:.2}% of hits begin in the initial or final 25% of target protein')
    
    return s, fpart
    

def gmgc_plot(s, fpart):
    import seaborn as sns
    from matplotlib import pyplot as plt
    
    fig, axarr = plt.subplot_mosaic([['a)', 'b)'], ['c)', 'd)'], ['e)', 'f)']], constrained_layout=True)
    print('-- plot identity vs. eval')
    sns.scatterplot(ax=axarr['a)'], data=s, x='pident', y='Log(E-value)', s=2, color='gray')
    axarr['a)'].set_xlabel('Identity (%)')
    axarr['a)'].set_ylabel('Log(E-value)')
    print('-- plot QxT len')
    sns.scatterplot(ax=axarr['b)'], data=s, x='qlen', y='tlen', s=5, color='gray')
    axarr['b)'].set_xlabel('AMP length (res.)')
    axarr['b)'].set_ylabel('Target length (res.)')
    print('-- plot Position of match starting')
    sns.kdeplot(ax=axarr['c)'], data=s, x='pct_start', hue='Identity', fill=True)
    axarr['c)'].set_xlabel('Match start (% of target length)')
    axarr['c)'].set_ylabel('Density (AU)')
    print('-- plot Qcov')
    sns.kdeplot(ax=axarr['d)'], data=s, x='qcov', hue='Identity', fill=True, legend=False)
    axarr['d)'].set_xlabel('Query coverage (%)')
    axarr['d)'].set_ylabel('Density (AU)')
    print('-- plot Tcov')
    sns.kdeplot(ax=axarr['e)'], data=s, x='tcov', hue='Identity', fill=True, legend=False)
    axarr['e)'].set_xlabel('Target coverage (%)')
    axarr['e)'].set_ylabel('Density (AU)')
    print('-- plot functions')
    sns.barplot(ax=axarr['f)'], data=fpart, x='functions', y='Log(counts)', color='gray')
    axarr['f)'].set_xticklabels(fpart.functions.tolist(),rotation = 30)
    axarr['f)'].set_xlabel('Functions')
    axarr['f)'].set_ylabel('Log(counts)')

    for label, ax in axarr.items():
        ax.set_title(label, fontfamily='Sans serif', loc='left', fontsize='large')

    plt.tight_layout()
    #plt.show()
    plt.savefig('out_test.svg')


def gmgc_plot2(s):
    import seaborn as sns
    from matplotlib import pyplot as plt
    
    fig, axarr = plt.subplots()
    sns.kdeplot(ax=axarr, data=s, x='pct_start', fill=False)
    axarr.set_xlabel('Match start (% of target length)')
    axarr.set_ylabel('Density (AU)')
    axarr.set_xlim(0, 100)
    axarr.set_yticks([])
    fig.savefig('homologs_matchstart.svg')


def classify(x):
    if 0 < x <= 25: return '0-25%'
    if 25 < x <= 50: return '25%-50%'
    if 50 < x <= 75: return '50%-75%'
    if 75 < x <= 100: return '75%-100%'


if __name__ == '__main__':
    main()
    
