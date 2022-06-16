def main():
    '''
    Create a heatmap of the shared environmental AMP contents
    '''
    import pandas as pd
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt
    from itertools import chain
    
    
    # creating sets of environments
    higher_level = {'sediment' : 'other',
        'bird gut' : 'other animal',
        'cat gut' : 'mammal gut',
        'insect associated' : 'other animal',
        'human urogenital tract' : 'other human',
        'dog gut' : 'mammal gut',
        'fermented food' : 'anthropogenic',
        'groundwater' : 'aquatic',
        'coral associated' : 'other animal',
        'rat gut' : 'mammal gut',
        'human associated' : 'other human',
        'cattle gut' : 'mammal gut',
        'deer gut' : 'mammal gut',
        'mouse gut' : 'mammal gut',
        'river associated' : 'aquatic',
        'primate gut' : 'mammal gut',
        'human respiratory tract' : 'other human',
        'cattle rumen' : 'other animal',
        'human saliva' : 'other human',
        'activated sludge' : 'anthropogenic',
        'lake associated' : 'aquatic',
        'wastewater' : 'anthropogenic',
        'chicken gut' : 'other animal',
        'air' : 'other',
        'human mouth' : 'other human',
        'plant associated' : 'soil/plant',
        'water associated' : 'aquatic',
        'pig gut' : 'mammal gut',
        'human skin' : 'other human',
        'marine' : 'aquatic',
        'soil' : 'soil/plant',
        'built environment' : 'anthropogenic',
        'human gut' : 'human gut',
        'anthropogenic': 'anthropogenic',
        'bear gut' : 'mammal gut',
        'rabbit gut': 'mammal gut',
        'algae associated': 'other',
        'crustacean gut': 'other animal',
        'cattle associated': 'other animal',
        'bird skin': 'other animal',
        'bee gut': 'other animal',
        'mussel associated': 'other animal',
        'fisher gut': 'mammal gut',
        'bat gut': 'mammal gut',
        'sponge associated': 'other animal',
        'human digestive tract': 'other human',
        'beatle gut': 'other animal',
        'dog associated': 'other animal',
        'insect gut': 'other animal',
        'extreme pH': 'other',
        'food': 'other',
        'guinea pig gut': 'mammal gut',
        'goat rumen': 'other animal',
        'mollusc associated': 'other animal',
        'goat gut': 'mammal gut',
        'horse gut': 'mammal gut',
        'wasp gut': 'other animal',
        'tunicate associated': 'other animal',
        'annelidae associated': 'other animal',
        'rodent gut': 'mammal gut',
        'ship worm associated': 'other animal',
        'coyote gut': 'mammal gut',
        'crustacean associated': 'other animal',
        'termite gut': 'other animal',
        'planarian associated': 'other animal',
        'thermal vent associated': 'other',
        'fish gut': 'other animal',
        'ice associated': 'other',
        'mock community': 'other',
        'mine': 'other',
        'pond associated': 'aquatic',
        'hot spring associated': 'other',
        }
    
    print('loading data')
    data = pd.read_table('data/gmsc_amp_genes_envohr_source.tsv.gz',
                         sep='\t',
                         header='infer')

    # filter duplicates
    data = data[data.is_metagenomic == True]
    data = data[['amp', 'general_envo_name']].drop_duplicates()
    data = data.groupby('general_envo_name')['amp'].apply(lambda x: set(x))

    # add environments with at least 100 peptides
    data = data.apply(lambda x: x if len(x) >= 100 else 'NA')
    data = data[~data.isin(['NA'])]

    # calculate overlap
    df = []
    for i in data.index:
        for j in data.index:
            n = len(data.loc[i, 'amp'].intersection(data.loc[j, 'amp']))
            df.append((i, j, n))
    
    df = pd.DataFrame(df, columns=['env1', 'env2', 'overlap'])
    df = df.pivot(index='env1', columns='env2', values='overlap')
    df.to_csv('lowlevel_habitat.tsv', sep='\t', header=True, index=True)
    
    # normalize
    df = df * 100 / df.max(axis=0)
    mask = np.zeros_like(df)
    mask[np.tril_indices_from(mask)] = True

    # plot heatmap
    sns.heatmap(df.astype('int'),
                annot=False,
                cmap="YlOrBr",
                square=True,
                xticklabels=True,
                yticklabels=True)
                
    plt.tight_layout()
    plt.xlabel(None)
    plt.ylabel(None)
    plt.savefig("low_level_habitats_ge100pep.svg")
    
    # convert to higher level
    data = data.reset_index()
    data['general_envo_name'] = data['general_envo_name'].apply(lambda x: higher_level.get(x, 'other'))
    data = data.groupby('general_envo_name')['amp'].apply(lambda x: set(chain.from_iterable(x)))

    # calculate overlap
    df = []
    for i in data.index:
        for j in data.index:
            n = data.loc[i].intersection(data.loc[j])
            n = len(n)
            df.append((i, j, n))

    df = pd.DataFrame(df, columns=['env1', 'env2', 'overlap'])
    df = df.pivot(index='env1', columns='env2', values='overlap')

    # normalize
    for c in df:
        df[c] = df[c]*100/df[c].max()
    
    df.round(2).to_csv('percent_bycolumn_overlap_higherenvo.tsv',
                       sep='\t',
                       header=True,
                       index=True)
                       
    # create mask of zeros
    mask = np.zeros_like(df)
    mask[np.tril_indices_from(mask)] = True

    # plot heatmap
    sns.heatmap(df.astype('int'),
                annot=False,
                cmap="YlOrBr",
                mask=mask, square=True)
    plt.tight_layout()
    plt.savefig('figure_1c_alt_heatmap_environments_overlap.svg')
    plt.close()


if __name__ == '__main__':
    main()
