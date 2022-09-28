def hrenvo():
    '''
    Creates the complete resource table used in the website
    as well as in many other analysis
    This table contains all gene info, taxonomy, specI and human
    readable environments
    '''
    import os
    import pandas as pd
    import numpy as np

    os.makedirs('outputs/',
                exist_ok=True)

    # load genes from ProGenomes2
    data = pd.read_table('data/AMPSphere_ProGenomes2.tsv.gz',
                         sep='\t',
                         header='infer')
    
    # load ncbi missing lineages
    missing_lineages = pd.read_table('data/ncbi_missing_lineages.txt')
    missing_lineages.rename({'taxname': 'source'}, axis=1, inplace=True)

    # import samples from progenomes
    refzero = pd.read_table('data/proGenomes2.1_specI_lineageNCBI.tab', sep='\t', header=None)
    refzero = refzero[[0, 7]]
    refzero = refzero.rename({0: 'genome',
                              7: 'source'},
                             axis=1)
                             
    ref = pd.read_table('data/progenomes_samples.tsv',
                        sep='\t',
                        header=None, 
                        names=['specI', 'genome'])

    ref = ref.merge(on='genome', right=refzero)
    
    source = []
    for x in ref['source'].tolist():
        if x == x:
            source.append(' '.join(x.split()[1:]))
        else:
            source.append(np.nan)
    
    ref['source'] = source    
    
    # merge samples, AMPs and genes from progenomes
    df1 = data.merge(on='genome', right=ref)
    
    # this step is needed because not all analyzed genomes
    # were in progenomes specI cluster tables
    df2 = data[~data.GMSC10.isin(df1.GMSC10)]
    df2['taxid'] = df2['genome'].apply(lambda x: x.split('.')[0]).astype('int')
    df2 = df2.merge(on='taxid', right=missing_lineages).drop('taxid', axis=1)
        
    # concatenating results and sorting it
    df = pd.concat([df1, df2])
    df = df.sort_values(by='AMP')

    # exporting progenomes complete table
    df.to_csv('outputs/pgenomes_AMP_specI.tsv.gz',
              sep='\t',
              header=True,
              index=None)

    # cleaning environment
    if len(df) == len(data): del data, ref, df1, df2, missing_lineages
    else:
        print('Failed in merging tables')
        break
        
    # load data from all metagenomic AMP genes
    df2 = pd.read_table('data/complete_amps_associated_taxonomy.tsv.gz',
                        sep='\t',
                        header='infer')

    # rename columns to compatibility
    df.rename({'GMSC10': 'gmsc',
               'AMP': 'amp',
               'genome': 'sample',
               'species': 'source'},
              axis=1,
              inplace=True)
              
    # shortening tables
    df = df[['gmsc', 'amp', 'sample', 'source', 'specI']]
    df2 = df2[['gmsc', 'amp', 'sample', 'source', 'specI']]

    # adding info about metagenome origins
    df['is_metagenomic'] = 'False'
    df = df.fillna('NA')
    df2['is_metagenomic'] = 'True'
    df2 = df2.fillna('NA')
    
    # concatenating results of AMPs originating from
    # progenomes and metagenomes
    gmsc_genes = pd.concat([df, df2])
    gmsc_genes = gmsc_genes.sort_values(by=['amp', 'gmsc'])
    gmsc_genes = gmsc_genes.reset_index(drop=True)
    gmsc_genes = gmsc_genes.replace('*', 'NA')

    # exporting table
    gmsc_genes.to_csv('outputs/complete_gmsc_pgenomes_metag.tsv.gz',
                      sep='\t',
                      header=True,
                      index=None)

    # cleaning environment
    del df, df2

    # loading environmental info
    metadata = pd.read_table('data/reduced_metadata.tsv',
                             sep='\t', 
                             header='infer')

    # renaming columns
    metadata.rename({'sample_accession':'sample'},
                    axis=1,
                    inplace=True)

    # merge metadata
    df = gmsc_genes.merge(on='sample', right=metadata)

    # selecting genes from progenomes,
    # which does not contain associated metadata
    df2 = gmsc_genes[~gmsc_genes.gmsc.isin(df.gmsc)]

    # concatenating genes from progenomes and metagenomes
    gdf = pd.concat([df, df2])
    gdf = gdf.sort_values(by=['amp', 'gmsc'])
    gdf = gdf.reset_index(drop=True)
    gdf = gdf.fillna('NA')
    gdf = gdf.replace('N.A.', 'NA')
    
    # export data
    gdf.to_csv('outputs/gmsc_amp_genes_envohr_source.tsv.gz',
               sep='\t',
               header=True,
               index=None)

