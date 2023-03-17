def filtermmseqs2_amp_contig():
    '''
    Filter the taxonomy annotated with mmseqs2 and GTDB
    The results of the mmseqs2 run are processed using
    the contig names to generate a smaller file
    that will be used later to annotate the AMPs
    '''
    import pandas as pd

    print('... loading metaG annotated AMPs')
    df = pd.read_table(f'data/AMPsphere_metaG_annotation.tsv.gz',
                       sep='\t', 
                       header='infer')

    print('... proceding to the filtering')
    fdf = pd.DataFrame()
    for j, dv in enumerate(pd.read_table('data/mmseqs2.lca_taxonomy.full.tsv.xz',
                                         sep='\t',
                                         header='infer',
                                         chunksize=10_000_000)):
        dv = df.merge(on=['sample', 'contig'], right=dv)
        fdf = pd.concat([fdf, dv], ignore_index=True)
        print(f'... processed chunk {j}')

    print('... recovering non-associated')
    nfdf = df[~df.gmsc.isin(fdf.gmsc)]
    fdf = pd.concat([fdf, nfdf])

    print('... sort values')
    fdf = fdf.sort_values(by=['accession', 'gmsc'])
    
    print('... data exporting')
    fdf.to_csv('data/gmsc_meta_taxo.tsv.gz',
               sep='\t',
               header=True,
               index=None)
        
