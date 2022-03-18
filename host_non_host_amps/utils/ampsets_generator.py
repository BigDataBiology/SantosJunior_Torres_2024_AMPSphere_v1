def select_amps(data):
    ## selecting amps
    import pandas as pd
    from .general_functions import aps
    print('... creating sets')
    mammal_amplist = data[data['general_envo_name'] == 'mammal']
    nonmammal_amplist = data[data['general_envo_name'] == 'non-mammal']
    plant_amplist = data[data['general_envo_name'] == 'plant']
    environmental_amplist = data[data['general_envo_name'] == 'environmental']

    mammal_amplist = set(mammal_amplist['amp'])
    nonmammal_amplist = set(nonmammal_amplist['amp'])
    plant_amplist = set(plant_amplist['amp'])
    environmental_amplist = set(environmental_amplist['amp'])

    print(f'''AMPs per set
    Mammal: {len(mammal_amplist)},
    Non-mammal: {len(nonmammal_amplist)},
    Plant: {len(plant_amplist)},
    Environmental: {len(environmental_amplist)}''')

    print('... outputting numbers')
    with open('analysis/amps_list_allsamples.txt', 'w') as ofile:
        ofile.write(f'Mammal: {mammal_amplist};\n')
        ofile.write(f'Non-mammal: {nonmammal_amplist};\n')
        ofile.write(f'Plant: {plant_amplist};\n')
        ofile.write(f'Environmental: {environmental_amplist};\n')
    
    print('... counting AMPs per sample per environment')
    d = pd.DataFrame()
    for name in ['mammal', 'non-mammal', 'plant', 'environmental']:
        d = pd.concat([d,
                       aps(data,
                           name)])
    
    d.reset_index(drop=True, inplace=True)

    return d


def get_sample_size(d):
    ## getting sample sizes
    import pandas as pd
    print('... loading sample info')    
    df_sample = pd.read_table('data/samples-min500k-assembly-prodigal-stats.tsv')
    print('... reducing sets')
    df_sample = df_sample[['sample_accession',
                           'assembly_total_length']]
    df_sample.rename({'sample_accession': 'sample'},
                     axis=1,
                     inplace=True)
    print('... merging elements')
    df = d.merge(on='sample',
                 right=df_sample)
    print('... normalization')
    df['normalized'] = df.AMPs * 1_000_000 / df.assembly_total_length
    print('... exporting')
    df.to_csv('analysis/all_amps_hostvsenv.tsv',
              sep='\t',
              header=True,
              index=None)
    
    return df


def ampsets(data):
    d = select_amps(data)
    df = get_sample_size(d)
    return df
    

