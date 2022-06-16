def main():
    '''
    Create a heatmap of overlapping AMP contents over 
    different human body sites
    '''
    import pandas as pd
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt

    print('Load data')
    
    data = pd.read_table('data/gmsc_amp_genes_envohr_source.tsv.gz',
                         sep='\t',
                         header='infer')

    # creating groups of AMPs per body site
    print('Filter data')
    skin = set(data[data['general_envo_name'] == 'human skin']['amp'])
    respiratory_tract = set(data[data['general_envo_name'] == 'human respiratory tract']['amp'])
    mouth = set(data[data['general_envo_name'].isin(['human mouth', 'human saliva'])]['amp'])
    digestive_tract = set(data[data['general_envo_name'] == 'human digestive tract']['amp'])
    gut = set(data[data['general_envo_name'] == 'human gut']['amp'])
    urogenital_tract = set(data[data['general_envo_name'] == 'human urogenital tract']['amp'])

    # create dictionary of names by position
    posix_dic = {0: 'skin',
                 1: 'respiratory_tract',
                 2: 'mouth',
                 3: 'digestive_tract',
                 4: 'gut',
                 5: 'urogenital_tract'}

    # calculating overlap
    matrix = [[],[],[],[],[],[]]
    setlists = [skin,
                respiratory_tract,
                mouth,
                digestive_tract,
                gut,
                urogenital_tract]

    print('Generating matrices')
    for n, i in enumerate(setlists):
        matrix[n].append(posix_dic[n])
        for m, j in enumerate(setlists):
            matrix[n].append(len(i.intersection(j)))

    # convert matrix into dataframe
    print('converting matrix')
    df = pd.DataFrame(matrix).set_index(0)
    df.columns = df.index
    df = df.astype('int')
    for c in df.columns:
        df[c] = df[c] * 100 / df[c].max()
    
    df.round(2).to_csv('percent_bycolumn_overlap_hbody.tsv',
                       sep='\t',
                       header=True,
                       index=True)
                       
    # creating mask of zeros
    mask = np.zeros_like(df)
    mask[np.tril_indices_from(mask)] = True

    # plot heatmap
    print('creating heatmaps')
    sns.heatmap(df.astype('int'),
                annot=False,
                cmap="YlOrBr",
                mask=mask,
                square=True)
                
    plt.tight_layout()
    plt.savefig('figure_1d_alt_heatmap_humanbodysites_overlap.svg')
    plt.close()
    

if __name__ == '__main__':
    main()
