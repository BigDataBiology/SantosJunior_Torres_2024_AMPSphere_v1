def family_size():
    '''
    Plot the size of families filtered by quality
    '''
    import pandas as pd
    import matplotlib.pyplot as plt
    from collections import Counter

    print('Load spheres')
    amp_fam = pd.read_table('data/SPHERE_v.2021-03.levels_assessment.tsv.gz', sep='\t', header='infer')
    amp_fam = amp_fam[['AMP accession', 'SPHERE_fam level III']]
    amp_fam = amp_fam.rename({'AMP accession': 'amp',
                              'SPHERE_fam level III': 'family'},
                             axis=1)
                             
    print('Load genes')
    # the file "taxonomy_annotation.tsv" is generated for Fig.1a
    data = pd.read_table('taxonomy_annotation.tsv')
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
    qual = pd.read_table('data/quality_families.txt', sep='\t', header='infer')
    
    # make sure families have a min size of 8 AMPs
    qual = qual[qual.total >= 8]
    
    # get the genera and families
    spec_fam = set(df.index)
    qual = set(qual.family)

    # plot pannel c (big plot)
    df.loc[qual.intersection(spec_fam)].hist(bins=100, grid=False)
    plt.ylabel('Quality-controlled families')
    plt.xlabel('Number of different genera')
    plt.savefig('qc_families_numberofgenera.svg')

    # plot pannel c (detailed plot)
    df.loc[qual.intersection(spec_fam)][df >= 20].hist(bins=100, grid=False)
    plt.xlabel('Number of different genera')
    plt.ylabel('Quality-controlled families')
    plt.savefig('qc_families_numberofgenera_ge20.svg')

