def plot_clusters(df, tag, z=None):
    import seaborn as sns
    import matplotlib.pyplot as plt
    if z == None: z = False
    metrics = ['euclidean', 'braycurtis', 'correlation']
    methods = ['single', 'complete', 'average', 'weighted']
    methods_euc = ['centroid', 'median', 'ward']
    for m in methods_euc:
        if z:
            sns.clustermap(df, cmap='YlOrBr', z_score=1, method=m, metric=metrics[0])
        else:
            sns.clustermap(df, cmap='YlOrBr', method=m, metric=metrics[0])
        plt.tight_layout()
        plt.savefig(f'{m}_euclidean_{tag}.svg')
        plt.close()
    for mt in methods:
        for me in metrics: 
            if z:
                sns.clustermap(df, cmap='YlOrBr', z_score=1, method=mt, metric=me)
            else:
                sns.clustermap(df, cmap='YlOrBr', method=mt, metric=me)
            plt.tight_layout()
            plt.savefig(f'{mt}_{me}_{tag}.svg')
            plt.close()
        
def load_merge():
    import pandas as pd
    from collections import Counter

    genes = pd.read_table('data/gmsc_amp_genes_envohr_source.tsv.gz')
    genes = genes[genes.is_metagenomic == True]
    genes = genes[['amp', 'sample', 'source', 'general_envo_name']]
    genes.general_envo_name = genes.general_envo_name.str.replace('human saliva', 'human mouth')
    
    # excluding environments with less than 100 AMPs
    keep = genes[['amp', 'general_envo_name']]
    keep = keep.drop_duplicates()
    keep = Counter(keep['general_envo_name'])
    keep = [k for k, v in keep.items() if v >= 100]
    genes = genes[genes.general_envo_name.isin(keep)]
    
    # counting amps per envo
    amp_per_envo = genes[['amp', 'general_envo_name']]
    amp_per_envo = amp_per_envo.drop_duplicates()
    amp_per_envo = Counter(amp_per_envo['general_envo_name'])
    amp_per_envo = dict(amp_per_envo)

    # counting amps per species
    genes = genes.fillna('*')
    genes.source = [x.split(' ')[0] for x in genes.source]
    amp_per_sp = genes[['amp', 'source']]
    amp_per_sp = amp_per_sp.drop_duplicates()
    amp_per_sp = Counter(amp_per_sp['source']).items()
    amp_per_sp = {k: v for k, v in amp_per_sp if v >= 100}
    
    # calculating samples by environment
    samples_by_envo = genes[['sample', 'general_envo_name']]
    samples_by_envo = samples_by_envo.drop_duplicates()
    samples_by_envo = samples_by_envo.groupby('general_envo_name')
    samples_by_envo = samples_by_envo.agg('size')
    samples_by_envo = samples_by_envo.to_dict()
    
    # calculating samples by species
    samples_by_sp = genes[['sample', 'source']]
    samples_by_sp = samples_by_sp.drop_duplicates()
    samples_by_sp = samples_by_sp.groupby('source')
    samples_by_sp = samples_by_sp.agg('size')
    samples_by_sp = samples_by_sp.to_dict()
    
    # merging motifs
    anno = pd.read_table('analysis/AMPSphere_v.2022-03.annotation.tsv')
    anno.rename({'id': 'amp'}, axis=1, inplace=True)
    genes = genes.merge(on='amp', right=anno)
    genes = genes[['amp', 'source',
                   'general_envo_name',
                   'motif_match']]
    
    # fixing missing fields and fixing source format
    genes = genes.fillna('w/o known motifs')
    
    return [genes, samples_by_sp, samples_by_envo,
            amp_per_sp, amp_per_envo]
           

def plot_species(genes, samples_by_sp, samples_by_envo,
                 amp_per_sp, amp_per_envo):
    import pandas as pd
    from itertools import chain
    from collections import Counter

    antimicrobial = {'γ-core motif', 'siderophore', 'ion-pair-pi',
                     'lfcinb', 'lipid binding motif', 'atcun',
                     'ngr motif', 'bact7', 'bfii',
                     'leulystrp', 'gala', 'bact5',
                     'gly-central–symmetrical', 'kld12', 'gly zipper',
                     'rana box'}

    ##############################################################
    # >>>  plotting top species <<<
    ##############################################################

    topsp = pd.DataFrame.from_dict(amp_per_sp, orient='index')
    topsp = topsp.sort_values(by=0)
    topsp = topsp.tail(30)
    topsp = topsp.drop(['Acidobacteriae', 'Burkholderiales', 'Actinobacteriota', 'Rhizobiales',
                        'Xanthobacteraceae', 'Bacteroidia', 'Oscillospiraceae',
                        'Bacteroidaceae', 'Oscillospirales', 'Firmicutes',
                        'Burkholderiaceae', 'Gammaproteobacteria', 'Bacteroidales',
                        'Actinomycetia', 'Lachnospiraceae', 'Alphaproteobacteria',
                        'Clostridia', 'Proteobacteria', 'unclassified',
                        'root', 'Bacteria', '*'], axis=0)
                        
    sp_motif = genes[genes.source.isin(topsp.index)]
    sp_motif = sp_motif.drop_duplicates()
    sp_motif = sp_motif[['amp', 'source', 'motif_match']]
    sp_motif.motif_match = [x.split('|') for x in sp_motif.motif_match]
      
    spmotif_df = pd.DataFrame()
    for record in sp_motif.groupby('source'):
        motif = pd.DataFrame.from_dict(Counter(chain.from_iterable(record[1].motif_match)), orient='index').T
        motif['source'] = record[0]
        spmotif_df = spmotif_df.append(motif)

    # eliminating motifs found in less than 100 AMPs
    spmotif_df = spmotif_df.fillna(0)
    spmotif_df.set_index('source', inplace=True)
    spmotif_df = spmotif_df.loc[:, spmotif_df.sum(axis=0) > 100]
    spmotif_df = spmotif_df.astype('int')
    
    # normalizing by AMPs and samples    
    newt = pd.DataFrame()
    for col in spmotif_df.T.columns:
        newt[col] = spmotif_df.T[col] / amp_per_sp[col]
        #newt[col] = newt[col] / samples_by_sp[col]  # normalize by number of samples
    
    # selecting just antimicrobial motifs                 
    sel = [x for x in antimicrobial if x in spmotif_df.columns]
    newt = newt.T
    newt = newt[sel]

    # plotting clustermaps with different algorithms
    # using normalized data of top_species represented
    # in AMPSphere
    newt.to_csv('top_sps_clustermap.tsv', sep='\t', header=True, index=True)
    plot_clusters(newt, 'top_species', True)


def plot_guts(genes, samples_by_sp, samples_by_envo,
              amp_per_sp, amp_per_envo):
    import pandas as pd
    from itertools import chain
    from collections import Counter

    antimicrobial = {'γ-core motif', 'siderophore', 'ion-pair-pi',
                     'lfcinb', 'lipid binding motif', 'atcun',
                     'ngr motif', 'bact7', 'bfii',
                     'leulystrp', 'gala', 'bact5',
                     'gly-central–symmetrical', 'kld12', 'gly zipper',
                     'rana box'}

    ##############################################################
    # >>>  plotting guts <<<
    ##############################################################
    gut_df = genes[genes.general_envo_name.str.contains('gut')]
    gut_df = gut_df[['amp', 'general_envo_name', 'motif_match']]
    gut_df = gut_df.drop_duplicates()
    gut_df.motif_match = [x.split('|') for x in gut_df.motif_match]
    gutdf = pd.DataFrame()
    for record in gut_df.groupby('general_envo_name'):
        motif = pd.DataFrame.from_dict(Counter(chain.from_iterable(record[1].motif_match)), orient='index').T
        motif['general_envo_name'] = record[0]
        gutdf = gutdf.append(motif)

    # fix missing fields
    gutdf = gutdf.fillna(0)
    gutdf.set_index('general_envo_name', inplace=True)
    gutdf = gutdf.loc[:, gutdf.sum(axis=0) > 100]
    gutdf = gutdf.astype('int')
    gutdf = gutdf[[x for x in antimicrobial if x in gutdf.columns]]

    # normalizing by AMPs and samples    
    newt = pd.DataFrame()
    for col in gutdf.T.columns:
        newt[col] = gutdf.T[col] / amp_per_envo[col]
        #newt[col] = newt[col] / samples_by_envo[col]  # normalize by number of samples
    
    # plotting clustermaps with different algorithms
    # using normalized data of gut environments represented
    # in AMPSphere
    newt = newt.T
    newt.to_csv('guts_clustermap.tsv', sep='\t', header=True, index=True)
    plot_clusters(newt, 'guts', True)
    

def plot_human(genes, samples_by_sp, samples_by_envo,
                 amp_per_sp, amp_per_envo):
    import pandas as pd
    from itertools import chain
    from collections import Counter

    antimicrobial = {'γ-core motif', 'siderophore', 'ion-pair-pi',
                     'lfcinb', 'lipid binding motif', 'atcun',
                     'ngr motif', 'bact7', 'bfii',
                     'leulystrp', 'gala', 'bact5',
                     'gly-central–symmetrical', 'kld12', 'gly zipper',
                     'rana box'}

    ##############################################################
    # >>>  plotting human sites <<<
    ##############################################################
    hdf = genes[genes.general_envo_name.str.contains('human')]
    hdf = hdf[['amp', 'general_envo_name', 'motif_match']]
    hdf = hdf.drop_duplicates()
    hdf.motif_match = [x.split('|') for x in hdf.motif_match]
    human = pd.DataFrame()
    for record in hdf.groupby('general_envo_name'):
        motif = pd.DataFrame.from_dict(Counter(chain.from_iterable(record[1].motif_match)), orient='index').T
        motif['general_envo_name'] = record[0]
        human = human.append(motif)

    # fix missing fields
    human = human.fillna(0)
    human.set_index('general_envo_name', inplace=True)
    human = human.loc[:, human.sum(axis=0) > 100]
    human = human.astype('int')
    human = human[[x for x in antimicrobial if x in human.columns]]

    # normalizing by AMPs and samples    
    newt = pd.DataFrame()
    for col in human.T.columns:
        newt[col] = human.T[col] / amp_per_envo[col]
        #newt[col] = newt[col] / samples_by_envo[col]  # normalize by number of samples
    
    # plotting clustermaps with different algorithms
    # using normalized data of gut environments represented
    # in AMPSphere
    newt = newt.T
    newt.to_csv('human_sites_clustermap.tsv', sep='\t', header=True, index=True)
    plot_clusters(newt, 'human_sites', True)


def plot_megaenvo(genes, samples_by_sp, samples_by_envo,
                  amp_per_sp, amp_per_envo):
    import pandas as pd
    from itertools import chain
    from collections import Counter
    from utils.environment_classification import environments_list
    
    antimicrobial = {'γ-core motif', 'siderophore', 'ion-pair-pi',
                     'lfcinb', 'lipid binding motif', 'atcun',
                     'ngr motif', 'bact7', 'bfii',
                     'leulystrp', 'gala', 'bact5',
                     'gly-central–symmetrical', 'kld12', 'gly zipper',
                     'rana box'}
    
    envdict = dict()
    for ev, vv in zip(['non_mammalian', 'mammalian_gut', 'mammalian_other', 
                       'built_env.', 'freshwater', 'extreme',
                       'wastewater', 'soil', 'plant associated',
                       'marine', 'other'], environments_list):
        for i in vv: envdict[i] = ev
    
    ##############################################################
    # >>>  plotting mega environments <<<
    ##############################################################
    
    motif_df = pd.DataFrame()
    mtable = genes[['amp', 'general_envo_name', 'motif_match']].drop_duplicates()
    mtable.motif_match = [x.split('|') for x in mtable.motif_match]
    for record in mtable.groupby('general_envo_name'):
        motif = pd.DataFrame.from_dict(Counter(chain.from_iterable(record[1].motif_match)), orient='index').T
        motif['envo'] = record[0]
        motif_df = motif_df.append(motif)
        
    motif_df = motif_df.fillna(0)
    motif_df = motif_df.set_index('envo')
    motif_df = motif_df.loc[:, motif_df.sum(axis=0) > 100]
    motif_df = motif_df.astype('int')
    motif_df = motif_df[[x for x in antimicrobial if x in motif_df.columns]]
    
    # fix 
    tochange = [k for k in envdict if 'human' in k]
    for k in tochange: envdict[k] = 'human other'
    envdict['human gut'] = 'human_gut'
    
    # normalizing by AMPs and samples    
    newt = pd.DataFrame()
    for col in motif_df.T.columns:
        newt[col] = motif_df.T[col] / amp_per_envo[col]
        #newt[col] = newt[col] / samples_by_envo[col]  # normalize by number of samples
    
    newt = newt.T
    newt.index = [ envdict[x] for x in newt.index ]
    newt = newt.reset_index().groupby('index').agg('mean')
    newt.to_csv('megaenvo_clustermap.tsv', sep='\t', header=True, index=True)
    plot_clusters(newt, 'mega_env', True)


def mplot():
    print('Loading items')
    genes, ssp, sse, aps, ape = load_merge()
    print('Plot top species')
    plot_species(genes, ssp, sse, aps, ape)
    print('Plot diff guts')
    plot_guts(genes, ssp, sse, aps, ape)
    print('Plot human body sites')
    plot_human(genes, ssp, sse, aps, ape)
    print('Plot meta environments')
    plot_megaenvo(genes, ssp, sse, aps, ape)   
    
