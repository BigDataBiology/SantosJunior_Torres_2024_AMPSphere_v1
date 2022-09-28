def get_info():
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
            'cattle rumen' : 'mammal gut',
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
            'bear gut' : 'mammal gut'}

    is_host_associated = {'human gut' : True,
            'soil/plant' : False,
            'aquatic' : False,
            'anthropogenic' : False,
            'other human' : True,
            'mammal gut' : True,
            'other animal' : True,
            'other' : False}

    return higher_level, is_host_associated


def envo(genes):
    higher_level, is_host_associated = get_info()
    envo = genes[['amp', 'general_envo_name']].drop_duplicates()
    envo.general_envo_name = envo.general_envo_name.apply(lambda g: higher_level.get(g, 'other'))
    envo['is_host_associated'] = envo.general_envo_name.apply(lambda x: is_host_associated.get(x, 'NA'))
    return envo


def motif_matrix(genes):
    from itertools import chain
    motifs = genes[['amp', 'motif_match']].drop_duplicates()
    motifs = motifs.set_index('amp')
    m = set(chain.from_iterable(motifs['motif_match'].apply(lambda x: x.split('|'))))
    for i in m:
        motifs[i] = motifs.motif_match.apply(lambda x: x.count(i))
    
    motifs = motifs.drop('motif_match', axis=1)
    motifs = motifs.reset_index()
    return motifs


def plot_umap_envo(envo, motifs):
    import umap
    import seaborn as sns
    import pandas
    import matplotlib.pyplot as plt
    reducer = umap.UMAP(random_state=42)
    env = envo.merge(on='amp', right=motifs).drop('amp', axis=1).drop_duplicates()
    env.reset_index(drop=True, inplace=True)
    reducer.fit(env[env.columns[2:]])
    embedding = reducer.transform(env[env.columns[2:]])
    principalDf = pd.DataFrame(data = embedding, columns = ['UMAP dim. 1', 'UMAP dim. 2'])
    finalDf = pd.concat([principalDf, env[env.columns[0:2]]], axis=1)
    sns.scatterplot(data=finalDf, x='UMAP dim. 1', y='UMAP dim. 2', hue='is_host_associated', s=5) 
    plt.savefig('host_associated.umap.svg')
    plt.close()
    sns.scatterplot(data=finalDf, x='UMAP dim. 1', y='UMAP dim. 2', hue='general_envo_name', s=5) 
    plt.savefig('higher_envo.umap.svg')
    plt.close()
    finalDf['is_human'] = finalDf.general_envo_name.str.contains('human')
    sns.scatterplot(data=finalDf, x='UMAP dim. 1', y='UMAP dim. 2', hue='is_human', s=5) 
    plt.savefig('human.umap.svg')
    plt.close()
    return env


def plot_umap_avg(genes, motifs):
    import umap
    import seaborn as sns
    import pandas
    import matplotlib.pyplot as plt
    reducer = umap.UMAP(random_state=42)
    df = genes.merge(on='amp', right=motifs)
    df = df.drop(['amp', 'source', 'motif_match'], axis=1)
    antimicrobial = ['γ-core motif', 'siderophore', 'ion-pair-pi',
                     'lfcinb', 'lipid binding motif', 'atcun',
                     'ngr motif', 'bact7', 'bfii',
                     'leulystrp', 'gala', 'bact5',
                     'gly-central–symmetrical', 'kld12', 'gly zipper',
                     'rana box']
    newdf = pd.DataFrame(columns=antimicrobial)
    for record in df.groupby('general_envo_name'):
        n = record[1][antimicrobial].sum(axis=0) / len(record[1])
        newdf.loc[record[0]] = n
    newdf = newdf.reset_index().rename({'index': 'env'}, axis=1)
    reducer.fit(newdf.drop('env', axis=1))
    embedding = reducer.transform(newdf.drop('env', axis=1))
    df = pd.DataFrame(embedding, columns=['umap component 1', 'umap component 2'])
    df = pd.concat([df, newdf['env']], axis=1)
    df['is_human'] = df.env.str.contains('human')
    higher_level, is_host_associated = get_info()
    df['higher'] = df.env.apply(lambda x: higher_level.get(x, 'other'))
    df['is_host_associated'] = df.higher.apply(lambda x: is_host_associated.get(x, 'other'))   
    sns.scatterplot(data=df, x='umap component 1', y='umap component 2', hue='is_host_associated')
    plt.savefig('avg_host_associated.umap.svg')
    plt.close()
    sns.scatterplot(data=df, x='umap component 1', y='umap component 2', hue='is_human')
    plt.savefig('avg_human.umap.svg')
    plt.close()
    sns.scatterplot(data=df, x='umap component 1', y='umap component 2', hue='higher')
    plt.savefig('avg_higher_environments.umap.svg')
    plt.close()


def main():
    from utils.motif_plot import load_merge
    genes, samples_by_sp, samples_by_envo, amp_per_sp, amp_per_envo = load_merge()
    motifs = motif_matrix(genes)
    env = envo(genes)
    plot_umap_envo(env, motifs)
    plot_umap_avg(genes, motifs)
    
    
if __name__ == '__main__':
    main()
    
