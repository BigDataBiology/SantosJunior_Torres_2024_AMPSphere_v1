def metadata():
    '''
    Creates a reduced and human-readable metadata table
    '''
    import pandas as pd
    metadata = pd.read_table('data/metadata.tsv.xz', sep='\t', header='infer', low_memory=False)
    metadata = metadata[['sample_accession',
                         'microontology',
                         'environment_material',
                         'geographic_location',
                         'latitude',
                         'longitude',
                         'host_scientific_name',
                         'host_tax_id']]

    envonames = pd.read_table('data/general_envo_names.tsv.xz', sep='\t', header='infer')
    envonames = envonames.fillna('N.A.')
    metadata = metadata.fillna('N.A.')
    metadata = metadata.merge(on=['microontology',
                                  'host_tax_id',
                                  'host_scientific_name'],
                              right=envonames)
    metadata = metadata[['sample_accession',
                         'geographic_location',
                         'latitude',
                         'longitude',
                         'general envo name',
                         'environment_material']]
    metadata.to_csv('data/reduced_metadata.tsv.xz',
                    sep='\t',
                    header=True,
                    index=None)
                    
