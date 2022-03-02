def heatmap_environments():
    '''
    Create a heatmap of the shared environmental AMP contents
    '''
    import pandas as pd
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    # creating sets of environments
    non_mammalian_host = ['bird gut', 'sponge associated', 'wasp gut',
                          'bee gut', 'planarian associated',  'crustacean gut',
                          'beatle gut', 'insect gut', 'mussel associated',
                          'termite gut', 'tunicate associated', 'bird skin',
                          'fish gut',  'ship worm associated', 'crustacean associated',
                          'annelidae associated', 'coral associated', 'mollusc associated',
                          'insect associated']
                          
    mammalian_host_guts = ['bear gut', 'coyote gut', 'pig gut',
                           'chicken gut', 'rodent gut', 'bat gut',
                           'fisher gut', 'rabbit gut', 'horse gut',
                           'rat gut', 'human gut', 'guinea pig gut',
                           'mouse gut', 'cattle gut', 'cat gut',
                           'deer gut', 'goat gut', 'primate gut',
                           'dog gut', 'dog associated']

    mammalian_host_other =  ['human respiratory tract', 'cattle associated',
                             'human mouth', 'human associated', 'human saliva',
                             'cattle rumen', 'human urogenital tract',
                             'goat rumen', 'human skin',
                             'human digestive tract']

    built_environment = ['built environment', 'anthropogenic']

    freshwater = ['ice associated', 'river associated',
                  'lake associated', 'groundwater',
                  'pond associated']
                  
    wastewater = ['wastewater', 'activated sludge']
    
    print('loading data')
    data = pd.read_table('data/gmsc_amp_genes_envohr_source.tsv',
                         sep='\t',
                         header='infer')

    # filtering data by environmental sets:
    non_mammal = set(data[data['general envo name'].isin(non_mammalian_host)]['amp'])
    mammal_gut = set(data[data['general envo name'].isin(mammalian_host_guts)]['amp'])
    mammal_others = set(data[data['general envo name'].isin(mammalian_host_other)]['amp'])
    built_environment = set(data[data['general envo name'].isin(built_environment)]['amp']) 
    freshwater = set(data[data['general envo name'].isin(freshwater)]['amp'])
    wastewater = set(data[data['general envo name'].isin(wastewater)]['amp']) 
    plant_associated = set(data[data['general envo name'] == 'plant associated']['amp']) 
    marine = set(data[data['general envo name'] == 'marine']['amp'])
    soil = set(data[data['general envo name'] == 'soil']['amp']) 

    # create dictionary of names by position
    posix_dic = {0: 'non_mammal',
                 1: 'mammal_gut',
                 2: 'mammal_others',
                 3: 'marine',
                 4: 'soil',
                 5: 'built_environment',
                 6: 'freshwater',
                 7: 'wastewater',
                 8: 'plant_associated'}

    # computing intersections over pairs of environments
    matrix = [[],[],[],[],[],[],[],[],[]]
    setlists = [non_mammal, mammal_gut, mammal_others, marine, soil, built_environment, freshwater, wastewater, plant_associated]
    for n, i in enumerate(setlists):
        matrix[n].append(posix_dic[n])
        for m, j in enumerate(setlists):
            matrix[n].append(len(i.intersection(j)))

    # convert overlap info into a dataframe
    df = pd.DataFrame(np.array(matrix)).set_index(0)
    df.columns = posix_dic.values()
    df = df.astype('int')

    # convert overlap into percent
    df = df * 100 / df.max()

    # create mask of zeros
    mask = np.zeros_like(df)
    mask[np.tril_indices_from(mask)] = True

    # plot heatmap
    sns.heatmap(df.astype('int'), annot=True, cmap="YlOrBr", mask=mask, square=True)
    plt.tight_layout()
    plt.savefig('figure_1c_heatmap_environments_overlap.svg')

