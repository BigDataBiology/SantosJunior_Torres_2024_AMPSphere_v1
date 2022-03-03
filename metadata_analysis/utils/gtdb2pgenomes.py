def conversion_from_GTDB_to_pGenomes():
    '''
    Create files to convert from GTDB taxonomy into 
    specI clusters
    '''
    import pandas as pd

    print('... load bacteria')

    # from GTDB
    data = pd.read_table('data/bac120_metadata_r202.tsv',
                         sep='\t',
                         header='infer',
                         low_memory=False)

    data = data[['accession','checkm_marker_lineage',
                 'gtdb_taxonomy','ncbi_assembly_name',
                 'ncbi_bioproject','ncbi_biosample',
                 'ncbi_organism_name','ncbi_taxid',
                 'ncbi_taxonomy']]
                 
    print('... load archaea')

    data2 = pd.read_table('data/ar122_metadata_r202.tsv',
                         sep='\t',
                         header='infer',
                         low_memory=False)

    data2 = data2[['accession','checkm_marker_lineage',
                 'gtdb_taxonomy','ncbi_assembly_name',
                 'ncbi_bioproject','ncbi_biosample',
                 'ncbi_organism_name','ncbi_taxid',
                 'ncbi_taxonomy']]
                 
    print('... load reference data')

    refdata = pd.read_table('data/progenomes_samples.tsv',
                            sep='\t',
                            header=None, 
                            names=['cluster', 'genome'])

    refdata['ncbi_taxid'] = [x[0] for x in refdata.genome.str.split('.')]
    refdata['ncbi_biosample'] = [x[1] for x in refdata.genome.str.split('.')]
    refdata.drop('genome', axis=1, inplace=True)

    df = pd.merge(on='ncbi_biosample',
                  right=refdata,
                  left=data).drop(['ncbi_tax_id'],
                                  axis=1)

    df2 = pd.merge(on='ncbi_biosample',
                   right=refdata,
                   left=data2).drop(['ncbi_tax_id'],
                                    axis=1)

    dfd = pd.concat([df, df2])

    print('... export 1')
    dfd.to_csv('data/gtdb_to_pgenomes.tsv',
               sep='\t',
               header=True,
               index=None)

    print('... cleaning')
    string_list = [ str(x).split(';') for x in df[['gtdb_taxonomy']].values]

    df['D'] = [ x[0] for x in string_list]
    df['P'] = [ x[1] for x in string_list]
    df['C'] = [ x[2] for x in string_list]
    df['O'] = [ x[3] for x in string_list]
    df['F'] = [ x[4] for x in string_list]
    df['G'] = [ x[5] for x in string_list]
    df['S'] = [ x[6] for x in string_list]

    df['D'] = [str(x).replace("['", '') for x in df.D.values ] 
    df['S'] = [str(x).replace("']", '') for x in df.S.values ] 

    df['D'] = [str(x).replace('d__', '') for x in df.D.tolist()]
    df['O'] = [str(x).replace('o__', '') for x in df.O.tolist()]
    df['F'] = [str(x).replace('f__', '') for x in df.F.tolist()]
    df['C'] = [str(x).replace('c__', '') for x in df.C.tolist()]
    df['P'] = [str(x).replace('p__', '') for x in df.P.tolist()]
    df['G'] = [str(x).replace('g__', '') for x in df.G.tolist()]
    df['S'] = [str(x).replace('s__', '') for x in df.S.tolist()]

    df = df.drop(['accession', 'checkm_marker_lineage',
                  'gtdb_taxonomy','ncbi_assembly_name',
                  'ncbi_bioproject', 'ncbi_organism_name',
                  'ncbi_taxid', 'ncbi_taxonomy'],
                 axis=1)

    print('... export 2')
    df.to_csv('data/conv_pgen_to_gtdb_list.tsv',
              sep='\t',
              header='infer',
              index=None)
              
