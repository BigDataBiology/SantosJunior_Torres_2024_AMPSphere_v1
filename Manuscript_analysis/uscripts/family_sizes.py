def main():
    '''
    Plot the size of families filtered by quality
    '''
    import pandas as pd
    import matplotlib.pyplot as plt
    from collections import Counter

    print('Load spheres')
    amp_fam = pd.read_table('data/SPHERE_v.2022-03.levels_assessment.tsv.gz', sep='\t', header='infer')
    amp_fam = amp_fam[['AMP accession', 'SPHERE_fam level III']]
    amp_fam = amp_fam.rename({'AMP accession': 'amp',
                              'SPHERE_fam level III': 'family'},
                             axis=1)
                             
    print('Load genes')
    # the file "taxonomy_annotation.tsv" is generated for Fig.1a
    data = pd.read_table('data/taxonomy_annotation.tsv.gz')
    data = data[['amp', 'fixed']].drop_duplicates()

    print('Convert amp into fams')
    data = data.merge(on='amp', right=amp_fam)
    data = data.sort_values('family')
    data = data.drop_duplicates()

    print('Calculating number of genera per family')
    gen_per_fam = data[['family', 'fixed']]
    gen_per_fam = gen_per_fam.drop_duplicates()
    gen_per_fam = gen_per_fam['family']
    gen_per_fam = Counter(gen_per_fam)    
    df = pd.DataFrame.from_dict(gen_per_fam,
                                orient='index',
                                columns=['number of genera'])

    df = df.sort_values(by='number of genera')

    print('Retrieving quality-controlled families')
    qual = pd.read_table('data/quality_families.txt.gz', sep='\t', header='infer')
    
    # make sure families have a min size of 8 AMPs
    f = (qual.experimental_evidence == True) | (qual.perc >= 75)
    qual = qual[f]
    qual = qual[(qual.total >= 8)]
    
    # plot pannel c (big plot)
    ndf = df.loc[df.index.isin(qual.family)]
    ndf = ndf.loc[(ndf <= 100)['number of genera']]
    ndf.hist(bins=100, grid=False)
    plt.ylabel('High quality families')
    plt.xlabel('Number of different genera')
    plt.xlim(0, 100)
    plt.savefig('qc_families_numberofgenera.svg')

    # plot pannel c (detailed plot)
    ndf = df.loc[df.index.isin(qual.family)]
    ndf = ndf.loc[(ndf > 100)['number of genera']]
    ndf.hist(bins=100, grid=False)
    plt.xlabel('Number of different genera')
    plt.ylabel('High quality families')
    plt.savefig('qc_families_numberofgenera_ge100.svg')

    print('Calculating top genera')
    d = data[data.family.isin(qual.family)]
    d = d[d.family.isin(df[df['number of genera'] == 1].index)]
    d = d[['family', 'fixed']]
    d = d.drop_duplicates()['fixed']
    d = pd.DataFrame.from_dict(Counter(d),
                               orient='index')

    d = d.sort_values(by=0)
    print(d.tail(10))


if __name__ == '__main__':
    main()
    
