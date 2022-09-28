def get_metadata():
    '''
    Retrieve metadata
    '''
    import pandas as pd
    df = pd.read_table('metadata/mag.taxonomy.gtdb.raw.tsv')
    df.rename({'user_genome': 'genome'}, axis=1, inplace=True)
    df.genome = [x.replace('.fa', '') for x in df.genome]
    df2 = pd.read_table('metadata/genome_data.MAGs.tsv')
    df2 = df2.merge(on='genome', right=df, how='outer')
    df2.to_csv('metadata/genomes_info_sebastian.tsv', sep='\t', header=True, index=None)
    return df2
       
       
def filter_df(metadata):
    '''
    filter for high-quality MAGs
    >input metadata: pandas dataframe
    '''
    c = (metadata.completeness >= 90)
    cx = (metadata.contamination <= 5)
    n = (metadata['n.contigs'] <= 500)
    n50 = (metadata['N50.contigs'] >= 10_000)
    h = (metadata['strain.heterogeneity'] <= 5)
    return (c & cx & n & n50 & h)


def macrel(genome, infolder, ofolder):
    import subprocess
    print(f'Macrel is processing {genome}')
    subprocess.call(['macrel', 
                     'contigs', 
                     '-f', f'{infolder}/{genome}.fa.gz',
                     '-o', f'{ofolder}/{genome}_macrel'])


def select_genomes():
    '''
    Select genomes by filtering the samples and retrieving clusters
    that appear in more than 1 type of subject
        
    ''' 
    import pandas as pd              
    metadata = pd.read_table(f'metadata/genomes_info_sebastian.tsv')
    metadata = metadata[filter_df(metadata)]
    xl_file = pd.read_excel(f'metadata/FMT_META.metadata.all.v01.xlsx',
                            sheet_name=False)
    xl_file = xl_file[['sample_alias',
                       'fmt.id',
                       'subject_type',
                       'timepoint.fmt',
                       'geographic_location',
                       'subject_disease_status',
                       'intervention',
                       'fmt_no',
                       'clinical_response']]
    xl_file.rename({'sample_alias': 'sample'},
                   axis=1,
                   inplace=True)
    df = metadata.merge(on='sample',
                        right=xl_file)
    df = df[['genome',
             'completeness',
             'contamination',
             'sample',
             'cluster',
             'subject_type',
             'timepoint.fmt',
             'clinical_response',
             'fmt.id']]
    df = df.sort_values(by=['subject_type',
                            'cluster',
                            'clinical_response',
                            'fmt.id'])
    df.to_csv(f'metadata/selected_genomes.tsv',
              sep='\t',
              header=True,
              index=None)
    return df
    
    
def count_cluster_per_condition(df):
    '''
    Count number of genomes happening in each condition
    per clusters
    '''    
    sel = df.pivot_table(index=['subject_type',
                                'cluster',
                                'clinical_response'],
                         aggfunc='size')
    sel = sel.reset_index()
    sel = sel.rename({0: 'counts'}, axis=1)
    sel.to_csv(f'metadata/counts_per_situation.tsv',
               sep='\t',
               header=True,
               index=None)


def soi(df):
    '''
    Export the list of species of interest
    '''
    from collections import Counter    
    with open(f'metadata/species_of_interest.txt', 'w') as ofile:
        for k, v in Counter(df.cluster).items():
            if v > 2:
                ofile.write(f'{k}\n')
    soi = [k for k, v in Counter(df.cluster).items() if v > 2]
    df = df[df.cluster.isin(soi)]
    return df


def find_amps(df):
    '''
    Run macrel over genomes of interest
    '''
    import os
    import glob
    from os.path import exists
    infile = []
    for f in glob.glob('macrel_results/*_macrel'):
        f = f.split('/')[-1].replace('_macrel', '')
        infile.append(f)
    for g in df.genome:
        if g in infile: pass
        else: macrel(g, 'data', 'macrel_results')


def clean():
    '''
    Eliminate unwanted genome files and leave
    just tar compressed file
    '''
    import os
    import glob
    for f in glob.glob('data/*.fa.gz'):
        os.remove(f)


def fmt_startup():
    print('Prepare metadata')
    get_metadata()
    print('Filtering genomes')
    df = select_genomes()
    print('Counting clusters per conditions')
    count_cluster_per_condition(df)
    print('Selecting species of interest')
    df = soi(df)
    print('Find AMPs')
    find_amps(df)
    print('Cleaning')
    clean()
