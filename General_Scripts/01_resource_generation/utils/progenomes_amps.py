def startup(data_folder, analysis_folder):
    '''
    Retrieve AMPs spotted in ProGenomes and 
    add it up to AMPSphere
    '''
    import pandas as pd
    from collections import Counter

    # loading data    
    data = pd.read_table(f'{data_folder}/GMSC10.proGenomesv2_05.AMPs.rename.tsv.gz',
                         sep='\t', header='infer', comment='#')

    ref = pd.read_table(f'{analysis_folder}/AMPsphere_GMSC_correspondence.tsv.gz',
                        sep='\t', header='infer', low_memory=False)

    ref = ref.rename({'sequence': 'Sequence'}, axis=1)
    
    specI = pd.read_table(f'{data_folder}/pgenomes_samples.tsv', 
                          sep='\t', header=None, names=['specI', 'genome'])

    listed_genomes = f'{data_folder}/combinedVsearch3.50p.20n.norm.opt.0.042.filteredBlackList.map.xz'
    listed_genomes = pd.read_table(listed_genomes,
                                   sep='\t', header=None,
                                   names=['specI', 'genome'])

    # adding progenomes info    
    data['genome'] = [ '.'.join(x[0:2]) for x in data['Access'].str.split('.')]
    data = data[data.genome.isin(listed_genomes.genome)]
    
    # selecting non_singleton AMPs
    non_singletons = set(ref['Sequence'])

    # cleaning and merging info    
    # here we add the pamps that match to metaG
    # and pamps that are present in more than 1 genome
    # but absent in metaG
    data = data[data['Sequence'].isin(non_singletons)] 
    df = ref.merge(on='Sequence', right=data[['Sequence', 'genome']])
    df = df.drop(['genes', 'n_of_genes', 'length'], axis=1)

    # adding specI
    sp = df.merge(on='genome', right=specI)
    nsp = data[~data.genome.isin(specI['genome'])][['Sequence', 'genome']]
    nsp['specI'] = '*'
    nsp = nsp.merge(on='Sequence', right=ref)
    nsp = nsp[['accession', 'Sequence', 'genome', 'specI']]

    df = pd.concat([sp, nsp])
    df = df.sort_values(by='accession')
    df = df.reset_index(drop=True)
    df.to_csv(f'{analysis_folder}/AMPsphere_proGenomes_correspondence.tsv.gz',
                sep='\t', header=True, index=None)
                

def link_taxa(data_folder, analysis_folder):
    '''
    Link NCBI taxonomy to the specI info
    '''
    import pandas as pd
    
    data = pd.read_table(f'{analysis_folder}/AMPsphere_proGenomes_correspondence.tsv.gz',
                         sep='\t', header='infer', low_memory=False)
                         
    ref = pd.read_table(f'{data_folder}/proGenomes2.1_specI_lineageNCBI.tab.xz',
                        sep='\t', header=None, names=['genome', 'kingdom',
                                                      'phylum', 'class',
                                                      'order', 'family',
                                                      'genus', 'species'])
    
    a = data.merge(on='genome', right=ref)
    b = data[~data.genome.isin(ref.genome)]
    df = pd.concat([a,b])
    df = df.fillna('*')
    df = df.sort_values(by='accession')
    df = df.reset_index(drop=True)
    df.to_csv(f'{analysis_folder}/AMPsphere_species_pGenomes.tsv.gz',
              sep='\t', header=True, index=None)

    
def ampsphere2progenomes():
    '''
    Link AMPs to their proGenomes taxonomy,
    genomes origins, and specI
    '''
    import os
    
    data_folder = 'data/'
    analysis_folder = 'analysis/'

    for d in [data_folder, analysis_folder]:
        os.makedirs(d, exist_ok=True)

    print('Select AMPs from ProGenomes')
    startup(data_folder, analysis_folder)
    
    print('Associating taxonomies')
    link_taxa(data_folder, analysis_folder)
    
