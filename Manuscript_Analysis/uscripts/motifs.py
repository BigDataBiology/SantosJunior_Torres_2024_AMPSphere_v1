def motif():
    motif = ['tle target', 'ion-pair-pi', 'vinculin',
             'gly-central–symmetrical', 'leulystrp',
             'cbp-a', 'dock_MAPK', 'guanine nucleotide exchange',
             'kld12', 'γ-core motif', 'golgi_lock',
             'arginine finger motif', 'cargo', 'walker',
             'bact7', 'cspa', 'kwarn', 'lfcinb', 'atcun',
             'pest-motif', 'siderophore', 'eno',
             'autoproteolytic cleavage motif', 'chemotaxis', 'spsa',
             'ngr motif', 'gly zipper', 'n-myristoylation/s-palmitoylation motif',
             'yada', 'ABC_motif', 'adherence', 'capsule', 'rgd',
             'lipid binding motif', 'gala', 'mtpm', 'rana box', 'type_iv',
             'bfii', 'bact5', 'caax box']
    return motif
    
    
def load_annotations():
    import pandas as pd
    from itertools import chain
    data = pd.read_table('data/AMPSphere_v.2022-03.annotation.tsv.gz',
                         sep='\t',
                         header='infer')
    data = data.fillna('wo_motif')
    data['motif_match'] = [x.split('|') for x in data.motif_match]
    data = data[['id', 'motif_match']]

    m = set(chain.from_iterable(data['motif_match']))
    for c in m:
        data[c] = [1 if c in x else 0 for x in data.motif_match]
    
    data = data.drop('motif_match', axis=1)    
    return data
    
    
def load_clusters():
    import pandas as pd
    data = pd.read_table('data/SPHERE_v.2022-03.levels_assessment.tsv.gz',
                         sep='\t', 
                         header='infer')
    data = data[['AMP accession', 'SPHERE_fam level III']]
    data = data.rename({'AMP accession': 'id',
                        'SPHERE_fam level III': 'family'},
                       axis=1)
    return data


def merge():
    anno = load_annotations()
    fam = load_clusters()
    return anno.merge(on='id', right=fam)


def famsize():
    import pandas as pd
    from collections import Counter
    fam = load_clusters()
    fam = fam.groupby('family').agg('size')
    fam = fam.reset_index()
    fam = fam.rename({'index': 'family',
                      0: 'seqs'},
                     axis=1)
    fam = fam.sort_values(by='family')
    return fam
    

def fixdfs():
    df_size = famsize()
    df = merge()
    f = df_size[df_size.seqs >= 8]['family']
    data = df[df.family.isin(f)]
    data = data.drop('id', axis=1)
    data = data.groupby('family')
    data = data.agg('sum')
    data.to_csv('motif_seqs_per_family_counts.tsv', 
                sep='\t',
                header=True,
                index=True)

    df_size = df_size[df_size.seqs >= 8]
    df_size = df_size.set_index('family')
    df_size = df_size['seqs']
    
    data = data.T / df_size
    data = data.T * 100
    data.to_csv('motif_seqs_per_family_norm.tsv', 
                sep='\t',
                header=True,
                index=True)

    return data

  
def countfams():
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt

    qual = pd.read_table('data/quality_families.txt.gz')
    qual = qual[qual.total >= 8]
    qual = qual[(qual.experimental_evidence == True) | (qual.perc >= 75)]
    qual = qual['family'].tolist()
    df = fixdfs()
    df = df.loc[qual]

    print('Processing motifs')
    tlist = []
    bdf = (df >= 75).T
    for c in bdf.columns:
        idx = ','.join(bdf[bdf[c]].index)
        tlist.append((c, idx))

    pd.DataFrame(tlist,
                 columns=['family',
                          'motifs']).to_csv('motifs_listed_per_family.tsv.gz', 
                                            sep='\t',
                                            header=True,
                                            index=None)
                 
    fampermotif = bdf.sum(axis=1).reset_index()
    fampermotif = fampermotif.sort_values(by=0)
    fampermotif = fampermotif[fampermotif[0] > 0]
    fampermotif = fampermotif.rename({'index': 'Motifs',
                                      0: 'High-quality families'},
                                     axis=1)

    fampermotif.to_csv('counts_motifs_per_family.tsv.gz', 
                       sep='\t', header=True, index=None)
                       
    fig, ax = plt.subplots()
    sns.barplot(ax=ax,
                data=fampermotif,
                y='Motifs',
                x='High-quality families',
                orient='h',
                color='darkblue')
    plt.tight_layout()
    fig.savefig('qc_fam_motifs.svg')


def main():
    countfams()


if __name__ == '__main__':
    main()
            
