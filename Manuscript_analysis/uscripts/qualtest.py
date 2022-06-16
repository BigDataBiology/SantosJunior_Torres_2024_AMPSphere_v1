def main():
    '''
    Assemble all results for separated quality tests
    '''
    import pandas as pd
    import matplotlib.pyplot as plt
    from collections import Counter

    print('... loading test results')
    df = pd.read_table('data/quality_assessment.tsv.gz', sep='\t', header='infer')
    df['Experimental evidence'] = 'Failed'
    df.loc[(df.metaproteomes == 'Passed') | (df.metatranscriptomes == 'Passed'), 'Experimental evidence'] = 'Passed'
    df.drop(['metaproteomes', 'metatranscriptomes'], axis=1, inplace=True)

    print('... counting amps per class')
    expevd = pd.DataFrame.from_dict(Counter(df['Experimental evidence']), orient='index').T
    antifam = pd.DataFrame.from_dict(Counter(df.Antifam), orient='index').T
    rna = pd.DataFrame.from_dict(Counter(df.RNAcode), orient='index').T
    terminal = pd.DataFrame.from_dict(Counter(df.Coordinates), orient='index').T
    
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
    plt.savefig('figure_S1a_amp_quality.svg')


if __name__ == '__main__':
    main()
    
