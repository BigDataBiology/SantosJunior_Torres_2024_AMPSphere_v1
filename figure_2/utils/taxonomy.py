def taxon_affiliation():
    '''
    Generate a table with the AMPs taxonomic affiliation
    which is useful to quick identify AMP origins througout 
    AMPSphere
    '''
    import pandas as pd
    import matplotlib.pyplot as plt
    from collections import Counter

    def selector(i: tuple):
        '''
        Classify the lowest taxonomical level associated to a given AMP
        '''
        if i[1]['level'].isin(['species']).sum() == 1: return f'{i[0]}\tspecies'
        if i[1]['level'].isin(['genus']).sum() == 1: return f'{i[0]}\tgenus'
        if i[1]['level'].isin(['family']).sum() == 1: return f'{i[0]}\tfamily'
        if i[1]['level'].isin(['order']).sum() == 1: return f'{i[0]}\torder'
        if i[1]['level'].isin(['class']).sum() == 1: return f'{i[0]}\tclass'
        if i[1]['level'].isin(['phylum']).sum() == 1: return f'{i[0]}\tphylum'
        if i[1]['level'].isin(['superkingdom']).sum() == 1: return f'{i[0]}\tsuperkingdom'
     
    # load table with the data
    data = pd.read_table('data/complete_amps_associated_taxonomy.tsv.gz', sep='\t', header='infer')

    # filter dataframe and drop duplicates
    dd = data[['amp', 'level', 'source', 'specI']].drop_duplicates()

    # getting annotated AMPs 
    annotated_amps = dd[~dd['level'].isna()]
    annotated_amps = annotated_amps['amp']
    annotated_amps = annotated_amps.drop_duplicates()
    print(f'# number of AMPs annotated at some level: {len(annotated_amps)}')

    # filter AMPs not associated to any origin
    ds = (~dd['level'].isna()) & (dd['level'] != 'no rank')
    ds = dd[ds]
    ds = ds[['amp', 'level']]
    ds = ds.drop_duplicates()
    ds = ds.sort_values(by=['amp', 'level'])
    ds = ds.reset_index(drop=True)

    # classify AMPs according lowest known levels   
    # and save results as a table
    with open('result_tax.txt', 'w') as ofile:
        ofile.write('amp\ttaxonomy\n')
        for i in ds.groupby('amp'):
            ofile.write(f'{selector(i)}\n')
            
    da = pd.read_table('result_tax.txt',
                       sep='\t',
                       header='infer')

    selected_group = da[da.taxonomy.isin(['species',
                                          'genus'])]

    selected_group = selected_group['amp']

    # convert dataframe of AMPs with annotated origins
    sp = data[(data['amp'].isin(selected_group)) & (data['level'].isin(['species', 'genus']))]
    sp = sp[['amp', 'taxid', 'level', 'source']]
    sp = sp.sort_values(by=['amp', 'taxid', 'source'])
    sp = sp.drop_duplicates()
    sp['fixed'] = [x.split(' ')[0] for x in sp.source.values]
    sp.to_csv('taxonomy_annotation.tsv', sep='\t', header=True, index=None)

    ## finding AMPs annotated to more than 1 genus:
    sp_next = sp[['amp', 'fixed']].drop_duplicates()

    notunigenus = Counter(sp_next.amp).values()
    notunigenus = Counter(notunigenus)    

    val = 0
    for k, v in notunigenus.items():
        if k > 1: val += v

    print(f'The number of AMPs actually associated to more than 1 genus is: {val}')

    ## finding 25 most common genera:
    ps = pd.DataFrame.from_dict(Counter(sp_next.fixed), orient='index')
    ps = ps.sort_values(0)
    ps = ps * 100 / 863498
    ps.tail(25).plot.bar(rot=30, legend=False)
    plt.ylabel('% of AMPSphere candidates')
    plt.xlabel('Top 25 genera in AMPSphere')
    plt.tight_layout()
    plt.savefig('top25genera.svg')

    ## finding 10 most common genera:
    ps.tail(10).plot.bar(rot=30, legend=False)
    plt.ylabel('% of AMPSphere candidates')
    plt.xlabel('Top 10 genera in AMPSphere')
    plt.tight_layout()
    plt.savefig('top10genera.svg')


def rank_taxa():
    '''
    Plot a bar chart with the lowest taxonomic level 
    identification per AMP
    '''
    import pandas as pd
    import matplotlib.pyplot as plt
    from collections import Counter

    data = pd.read_table('result_tax.txt')

    # creating 
    d = dict(Counter(data.taxonomy.str.capitalize()))
    none = 863498 - sum(d.values())  # 863,498 is the total number of AMPs in AMPSphere
    d['None'] = none

    df = pd.DataFrame.from_dict(d,
                                orient='index',
                                columns=['AMP candidates']).loc[['None',
                                                                 'Superkingdom',
                                                                 'Phylum',
                                                                 'Class',
                                                                 'Order',
                                                                 'Family',
                                                                 'Genus',
                                                                 'Species']]

    df = df * 100 / df.sum()

    df.T.plot.bar(stacked=True, cmap='Dark2')
    plt.tight_layout()
    plt.ylabel('AMPSphere candidates - %')
    plt.savefig('taxonrank.svg')


def taxon_analysis():
    print('Taxon affiliation')
    taxon_affiliation()
    print('Ranking the taxa by AMP numbers')
    rank_taxa()
    
