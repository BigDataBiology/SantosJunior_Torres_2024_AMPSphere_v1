def addspecI():
    '''
    Add specI cluster info in the complete table of AMP, genes, and taxonomy
    '''
    import gzip
    import pandas as pd
    from collections import Counter

    # define taxonomy levels
    tab_levels = {'genus': 'G',
                  'order': 'O',
                  'species': 'S',
                  'superkingdom': 'D',
                  'class': 'C',
                  'family': 'F',
                  'phylum': 'P'}

    print('... loading amps')
    df = pd.read_table('data/gmsc_meta_taxo.tsv.gz',
                       sep='\t',
                       header='infer')

    print('... loading reference')
    reftab = pd.read_table('data/conv_pgen_to_gtdb_list.tsv',
                           sep='\t',
                           header='infer')
                           
    reftab = reftab.rename({'ncbi_biosample': 'sample'}, axis=1)
    
    for i in reftab.columns:
        if i != 'cluster':
            reftab[i] = reftab[i].str.replace('_', ' ')
    
    print('... processing pairs')
    # perfect sample match
    directly_assoc = df.merge(on='sample', right=reftab[['sample', 'cluster']])
    
    # impossible match due to no rank or level n.a.n.
    non_assoc = df[(df.level.isna()) | (df.level == 'no rank')]
    non_assoc['cluster'] = '*'
    
    # reducing the data already done
    df = df[~df.gmsc.isin(directly_assoc.gmsc)]
    df = df[~df.gmsc.isin(non_assoc.gmsc)]
    
    # associating specI by name and tax level
    # creating the sets per tax level of singled specI names
    selected_taxa = dict()
    for i in tab_levels.values():
        small_list = reftab[[i, 'cluster']]
        small_list = small_list.drop_duplicates()
        small_list = Counter(small_list[i]).items()
        selected_taxa[i] = [k for k,v in small_list if v == 1]

    newdf = pd.DataFrame()
    for i in df.groupby(['level', 'name']):
        print(f'{i[0]}')
        # get set by tax level
        v = tab_levels[i[0][0]]
        if i[0][1] in selected_taxa[v]:
            # get set by tax name
            g = (reftab[v] == i[0][1])
            # retrieve cluster
            g = reftab[g]['cluster']
            g = g.drop_duplicates()
            g = g.tolist()
            # eliminate multi-cluster results
            i[1]['cluster'] = g[0]
            # add results to previous frames
            newdf = pd.concat([newdf, i[1]])
          

    df = df[~df.gmsc.isin(newdf.gmsc)]
    df['cluster'] = '*'
    
    # uniting all results
    df = pd.concat([directly_assoc,
                    non_assoc,
                    newdf,
                    df])
                    
    df = df.rename({'cluster': 'specI',
                    'accession': 'amp',
                    'gene': 'gmsc',
                    'name': 'source'},
                   axis=1)
    
    df = df.sort_values(by=['amp',
                            'gmsc',
                            'level',
                            'source'])
                                
    print('... saving table')
    df.to_csv('data/complete_amps_associated_taxonomy.tsv.gz',
              sep='\t', 
              header=True,
              index=None)     
              
