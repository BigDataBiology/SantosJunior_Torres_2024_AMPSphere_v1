def metadata():
    '''
    Creates a reduced and human-readable metadata table
    '''
    import pandas as pd

    metadata = pd.read_table('data/metadata.tsv', sep='\t', header='infer', low_memory=False)

    metadata = metadata[['sample_accession',
                         'microontology',
                         'environment_material',
                         'geographic_location',
                         'latitude',
                         'longitude',
                         'host_scientific_name',
                         'host_tax_id']]

    metadata.microontology = [x.replace('-', ' ').replace('_', ' ') for x in metadata.microontology.tolist()]

    envonames = pd.read_table('data/general_envo_names.tsv', sep='\t', header='infer')

    envonames = envonames[['microontology',
                           'host scientific name',
                           'host tax id',
                           'general envo name']]

    envonames.rename({'host scientific name': 'host_scientific_name',
                      'host tax id': 'host_tax_id'},
                     axis=1,
                     inplace=True)

    metadata = metadata.merge(on=['microontology',
                                  'host_tax_id', 'host_scientific_name'],
                              right=envonames)

    metadata = metadata[['sample_accession',
                         'geographic_location',
                         'latitude',
                         'longitude',
                         'general envo name',
                         'environment_material']]

    metadata.to_csv('data/reduced_metadata.tsv',
                    sep='\t',
                    header=True,
                    index=None)
                    
