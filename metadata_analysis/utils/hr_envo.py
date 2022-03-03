def hrenvo():
    '''
    Creates the complete resource table used in the website
    as well as in many other analysis
    This table contains all gene info, taxonomy, specI and human
    readable environments
    '''
    import pandas as pd
    
    # imports genes from progenomes
    data = pd.read_table('data/AMPSphere_ProGenomes2.tsv.gz',
                         sep='\t',
                         header='infer')

    # import samples from progenomes
    ref = pd.read_table('data/pgenomes_samples.tsv',
                        sep='\t',
                        header=None, 
                        names=['specI', 'genome'])

    # merge samples, AMPs and genes from progenomes
    df1 = data.merge(on='genome', right=ref)
    
    # this step is needed because not all analyzed genomes
    # were in progenomes specI cluster tables
    df2 = data[~data.AMP.isin(df1.AMP)]
    
    # concatenating results and sorting it
    df = pd.concat([df1,df2])
    df = df.sort_values(by='AMP')

    # exporting progenomes complete table
    df.to_csv('pgenomes_AMP_specI.tsv',
              sep='\t',
              header=True,
              index=None)

    # cleaning environment
    del data, ref, df1, df2

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
              
    df2.rename({'name': 'source'},
               axis=1,
               inplace=True)

    # shortening tables
    df = df[['gmsc', 'amp', 'sample', 'source', 'specI']]
    df2 = df2[['gmsc', 'amp', 'sample', 'source', 'specI']]

    # adding info about metagenome origins
    df['is_metagenomic'] = False
    df2['is_metagenomic'] = True

    # concatenating results of AMPs originating from
    # progenomes and metagenomes
    gmsc_genes = pd.concat([df, df2])
    gmsc_genes = gmsc_genes.fillna('*')
    gmsc_genes = gmsc_genes.sort_values(by='gmsc')
    gmsc_genes = gmsc_genes.reset_index(drop=True)

    # exporting table
    gmsc_genes.to_csv('complete_gmsc_pgenomes_metag.tsv',
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
    gdf = gdf.sort_values(by='gmsc')
    gdf = gdf.reset_index(drop=True)

    # export data
    gdf.to_csv('gmsc_amp_genes_envohr_source.tsv',
               sep='\t',
               header=True,
               index=None)

