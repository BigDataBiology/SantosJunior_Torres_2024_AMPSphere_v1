def heatmap_mammal_guts():
    '''
    Caculate heatmaps of AMP overlap contents in 
    different mammalian guts
    '''
    import pandas as pd
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt
    

    print('Load data')
    data = pd.read_table('data/gmsc_amp_genes_envohr_source.tsv.gz',
                         sep='\t',
                         header='infer')

    print('# getting AMP lists for each host')
    human_gut = set(data[data['general envo name'].isin(['human gut'])]['amp'].tolist())
    pig_gut = set(data[data['general envo name'].isin(['pig gut'])]['amp'].tolist())
    chicken_gut = set(data[data['general envo name'].isin(['chicken gut'])]['amp'].tolist())
    mouse_gut = set(data[data['general envo name'].isin(['mouse gut'])]['amp'].tolist())
    dog_gut = set(data[data['general envo name'].isin(['dog gut'])]['amp'].tolist())
    cat_gut = set(data[data['general envo name'].isin(['cat gut'])]['amp'].tolist())
    bovine_gut = set(data[data['general envo name'].isin(['cattle gut'])]['amp'].tolist())

    # stating hosts
    posix_dic = {0: 'human_gut',
                 1: 'pig_gut',
                 2: 'chicken_gut',
                 3: 'mouse_gut',
                 4: 'cat_gut',
                 5: 'dog_gut',
                 6: 'bovine_gut'}

    # calculating overlaps
    matrix = [[],[],[],[],[],[],[]]
    setlists = [human_gut, pig_gut, chicken_gut, mouse_gut, cat_gut, dog_gut, bovine_gut]

    for n, i in enumerate(setlists):
        matrix[n].append(posix_dic[n])
        for m, j in enumerate(setlists):
            matrix[n].append(len(i.intersection(j)))

    # converting overlap info into a dataframe
    df = pd.DataFrame(np.array(matrix)).set_index(0)
    df.columns = posix_dic.values()
    df = df.astype('int')

    # converting overlap df into percent
    df = df * 100 / df.max()
    df = df.astype('int')

    # creat mask of zeros
    mask = np.zeros_like(df)
    mask[np.tril_indices_from(mask)] = True

    # plot heatmap
    sns.heatmap(df, annot=True, cmap='YlOrBr', mask=mask, square=True)
    plt.tight_layout()
    plt.savefig('figure_S1d_heatmap_mammalguts_overlap.svg')
    
