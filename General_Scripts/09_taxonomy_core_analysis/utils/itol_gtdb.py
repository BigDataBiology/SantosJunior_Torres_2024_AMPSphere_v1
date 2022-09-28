def itol_gtdb_prep():
    import pandas as pd
    from Bio import Phylo
    from re import sub
    from tqdm import tqdm
    from matplotlib import pyplot as plt
    from matplotlib.colors import to_hex
    
    # load lineages from GTDB
    data = pd.read_table('data/bac120_taxonomy_r202.tsv', names=['genome', 'lineage'])
    data['genus'] = data.lineage.apply(lambda x: sub('^.*;g__', '', x))
    data['genus'] = data.genus.apply(lambda x: sub(';.*', '', x))
    
    # load tree
    gtdb_tree = Phylo.read('data/bac120_r202.tre', 'newick')
    data = data[data.genome.isin([i.name for i in gtdb_tree.get_terminals()])]
    
    # load amps per assembled gigabase pairs
    ampspergbp = pd.read_table('analysis/normalized_amps_per_bp_per_taxon.tsv.xz')
    ampspergbp.rename({'fixed': 'genus'}, axis=1, inplace=True)
    ampspergbp = ampspergbp[ampspergbp.VAR_pct <= 10]
    
    # merge information
    data = data.merge(on='genus', right=ampspergbp)
    
    # keep only 1 genome per taxon
    ndata = data.groupby('genus').apply(lambda x: x.head(1))
    ndata = ndata.reset_index(drop=True)
    
    # create a hash table for genomes and genera 
    taxkeys = ndata[['genome', 'genus']]
    taxkeys = taxkeys.set_index('genome')
    taxkeys = taxkeys.to_dict()['genus']
    
    # eliminate genomes without AMP density
    # or redundant
    for i in tqdm(gtdb_tree.get_terminals()):
        if i.name in taxkeys: i.name = taxkeys[i.name]
        else: gtdb_tree.prune(i)
    
    # excluding discrepant value (5 times the second most abundant)
    # gtdb_tree.prune('Shimwellia')
    
    # export tree file for iTOL
    Phylo.write(gtdb_tree,
                'analysis/out.tre',
                'newick')
    
    # export annotation of AMP density
    ndata[['genus',
           'amps_per_Gbp']].to_csv('analysis/out.anno.txt',
                                   sep='\t',
                                   header=None,
                                   index=None)
    
    # reduce tree
    # tree is too big with many small phylum
    # reducing tree
    ndata['phylum'] = ndata['lineage'].apply(lambda x: x.split(';')[1][3:])
    k = ndata.groupby('phylum').agg('size')
    k = k[k >= 15].index
    ndata = ndata[ndata.phylum.isin(k)]
    ksp=set(ndata.genus)
    
    # edit tree 
    for i in gtdb_tree.get_terminals():
        if i.name not in ksp: gtdb_tree.prune(i.name)
    
    Phylo.write(gtdb_tree,
                'analysis/red_out.tre',
                'newick')
    
    # export annotation of phylum
    cmap = plt.cm.get_cmap('tab20')
    phylum = dict()
    for idx, i in enumerate(set(ndata.phylum)):
        phylum[i] = to_hex(cmap(idx))
    
    # color keys
    # {'Acidobacteriota': '#1f77b4', 'Actinobacteriota': '#aec7e8',
    #  'Bacteroidota': '#ff7f0e', 'Firmicutes': '#ffbb78',
    #  'Firmicutes_A': '#2ca02c', 'Proteobacteria': '#98df8a',
    #  'Verrucomicrobiota': '#d62728'}
    
    ndata['color'] = ndata.phylum.apply(lambda x: phylum[x])
    ndata[['genus',
           'color',
           'phylum']].to_csv('analysis/out.anno_ph.txt',
                             sep='\t',
                             header=None,
                             index=None)
    
