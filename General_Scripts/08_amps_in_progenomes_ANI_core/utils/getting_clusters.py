def get_clusters(data, filter):
    import networkx as nx
    nodes = []
    for _, s1, s2, ani, _, _ in data.loc[filter].itertuples():
        nodes.append((s1, s2))
    G = nx.Graph()
    G.add_edges_from(nodes)
    return list(nx.connected_components(G))


def format_labels(labels):
    newdict = dict()
    for idx, i in enumerate(labels):
        for j in i: newdict[j] = idx
    return newdict


def get_clones_strains(infile, cutoffstrain, cutoffclone):
    import pandas as pd
    data = pd.read_table(infile,
                         header=None,
                         names=['s1', 's2', 'ani', 'match', 'total'])
    data.loc[data.s1 == data.s2, 'ani'] = 100
    data.s1 = data.s1.str.replace('.fna.gz', '', regex=False)
    data.s2 = data.s2.str.replace('.fna.gz', '', regex=False)
    c = get_clusters(data,
                     (data.ani >= cutoffclone))
    st = get_clusters(data,
                      (data.ani >= cutoffstrain))
    c = pd.DataFrame.from_dict(format_labels(c), orient='index')
    c = c.rename({0: 'clone'}, axis=1)
    st = pd.DataFrame.from_dict(format_labels(st), orient='index')
    st = st.rename({0: 'strain'}, axis=1)
    return pd.concat([c, st], axis=1)


def load_ampsphere():
    import pandas as pd
    df = pd.read_table('data/complete_gmsc_pgenomes_metag.tsv.gz')
    df = df[df.is_metagenomic==False]
    df['sample'] = df['sample'].apply(lambda x: '.'.join(x.split('.')[1:]))
    df = df.groupby('sample').apply(lambda x: set(x.amp))
    return df


def head_row(w):
   w = w.reset_index()
   w = w['index'].sample(1)
   return w


def testamps(amps, s1, s2):
    if (s1 in amps.index): as1 = amps.loc[s1]
    else: as1 = set()    
    if (s2 in amps.index): as2 = amps.loc[s2] 
    else: as2 = set()
    overlap = as1.intersection(as2)
    total = as1.union(as2)
    no, ntotal = len(overlap), len(total)
    if ntotal == 0: npct=0
    else: npct = no*100/ntotal
    return (no, ntotal, npct)


def test_overlap(infile, amps, cutoffstrain, cutoffclone):
    import pandas as pd
    from itertools import combinations
    name = infile.split('/')[-1]
    name = name.replace('_ANI.tsv', '')
    df = get_clones_strains(infile, cutoffstrain, cutoffclone)
    df.to_csv(f'amp_results/{name}_strain_clones.tsv',
              sep='\t',
              header=True,
              index=True)    
    strain_sample = df.groupby('clone').apply(lambda w: head_row(w))
    strain_sample = df.loc[strain_sample[0].tolist(),
                           'strain']
    strain_sample = strain_sample.reset_index()
    strain_sample = strain_sample.rename({'index': 'genome'}, axis=1)
    strain_sample = strain_sample.groupby('strain')
    strain_sample = strain_sample.apply(lambda x: x.genome.tolist())
    strain_sample = strain_sample.reset_index().rename({0: 'genomes'}, axis=1)
    strain_sample['L'] = strain_sample.genomes.apply(lambda x: len(x))
    overl=[]
    if len(strain_sample[strain_sample.L > 1]) > 0:
        print('Computing within-strain overlaps')
        for st in strain_sample.loc[strain_sample.L > 1, 'strain']:
            genomes = strain_sample.loc[strain_sample.strain == st, 'genomes'].tolist()[0]
            for g1, g2 in combinations(genomes, 2):
                no, ntotal, npct = testamps(amps, g1, g2)    
                overl.append((name, g1, g2, st, st, no, ntotal, npct))      
    print('Computing cross-strain overlaps')
    for st1, st2 in combinations(strain_sample.strain, 2):
        genomes1 = strain_sample.loc[strain_sample.strain == st1, 'genomes'].tolist()[0]
        genomes2 = strain_sample.loc[strain_sample.strain == st2, 'genomes'].tolist()[0]
        for g1 in genomes1:
            for g2 in genomes2:
                no, ntotal, npct = testamps(amps, g1, g2)    
                overl.append((name, g1, g2, st1, st2, no, ntotal, npct))
    overl = pd.DataFrame(overl, columns=['species', 'genome1', 'genome2',
                                         'strain1', 'strain2', 'shared_amps',
                                         'total_nr_amps', 'percent_shared'])
    overl.to_csv(f'amp_results/{name}_amp_overlaps.tsv',
                 sep='\t',
                 header=True,
                 index=None)       


def analyze_pairs():
    import pandas as pd
    from glob import glob
    res = []
    for infile in glob('amp_results/*_amp_overlaps.tsv'):
        name = infile.split('/')[-1]
        name = name.replace('_amp_overlaps.tsv', '')
        df = pd.read_table(infile)
        ss_share = len(df[(df.strain1 == df.strain2) & (df.shared_amps != 0)])
        ss_nshare =  len(df[(df.strain1 == df.strain2) & (df.shared_amps == 0)])
        ds_share = len(df[(df.strain1 != df.strain2) & (df.shared_amps != 0)])
        ds_nshare =  len(df[(df.strain1 != df.strain2) & (df.shared_amps == 0)])
        res.append((name, ss_share, ss_nshare,
                    ds_share, ds_nshare))
    res = pd.DataFrame(res, columns=['species', 'pairs_of_same_strain_sharing_AMPs',
                                     'pairs_of_same_strain_NOT_sharing_AMPs', 
                                     'pairs_of_dif_strain_sharing_AMPs',
                                     'pairs_of_dif_strain_NOT_sharing_AMPs'])
    return res
                    

def test_ani_values():
    import pandas as pd
    from glob import glob
    from scipy.stats import percentileofscore as pos
    reslist = []
    for infile in glob('analysis/*.tsv'):
        name = infile.split('/')[-1].replace('_ANI.tsv', '')
        df = pd.read_table(infile, names=['s1', 's2', 'ani', 'm', 't'])
        w = pos(df.ani, 99.5)
        x = pos(df.ani, 99.99)
        mean = df.ani.mean()
        std = df.ani.std()
        q25, q50, q75, q90, q95, q99 = df.ani.quantile([.25,.5,.75,.90,.95,.99])
        reslist.append((name, w, x, mean, std, q25, q50, q75, q90, q95, q99))
    reslist = pd.DataFrame(reslist, columns=['species', 'percentileofANI99.5',
                                             'percentileofANI99.99', 'mean',
                                             'std', 'q25', 'q50', 'q75', 'q90',
                                             'q95', 'q99'])
    return reslist

                    
def main(standard=None):
    import os
    from glob import glob
    from scipy.stats import fisher_exact
    if standard == None:
        standard = True
    print('Making output dir')
    os.makedirs('amp_results', exist_ok=True)
    print('Loading AMPSphere')
    amps = load_ampsphere()
    print('Determining cutoffs')
    vals = test_ani_values()
    vals.to_csv('amp_results/ANI_analysis_by_species.tsv',
                sep='\t', header=True, index=None)
    vals = vals.set_index('species')
    vals = vals[['q75', 'q90']]
    vals = vals.to_dict()            
    for infile in glob('analysis/*.tsv'):
        print(f'Analyzing {infile}')
        if standard:
            test_overlap(infile, amps, 99.5, 99.99)
        else:
            name = infile.split("/")[-1].replace('_ANI.tsv', '')
            cutoffstrain = vals['q75'][name]
            cutoffclone = vals['q90'][name]
            test_overlap(infile, amps, cutoffstrain, cutoffclone)
    print('Testing probabilities')
    df = analyze_pairs()
    df.to_csv('pairs_same_dif_strains.tsv',
              sep='\t',
              header=True,
              index=None)
    df = df.set_index('species')
    df = df[df.sum(axis=1) != 0]
    df = df.sum(axis=0)
    p1 = df.pairs_of_same_strain_sharing_AMPs + df.pairs_of_same_strain_NOT_sharing_AMPs
    p1 = df.pairs_of_same_strain_sharing_AMPs * 100 / p1
    p2 = df.pairs_of_dif_strain_sharing_AMPs + df.pairs_of_dif_strain_NOT_sharing_AMPs
    p2 = df.pairs_of_dif_strain_sharing_AMPs * 100 / p2
    print(f'''The proportion of getting a pair of genomes sharing AMPs
    from the same strain is {p1:.2f} and from different strains
    is {p2:.2f} controling for the species during the comparison''')
    table = [df.tolist()[0:2],
             df.tolist()[2:4]]
    odds, p = fisher_exact(table, alternative='two-sided')
    print(f'''The odds are {odds:.1f}-fold (P={p:.2E}) higher of having a pair of
    genomes belonging to the same strain sharing AMPs than a pair of genomes
    from different strains in the same species sharing AMPs''')         


if __name__ == '__main__':
    main()
    

if __name__ == '__main__':
    main()
