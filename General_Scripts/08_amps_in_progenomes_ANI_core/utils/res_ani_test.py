def checkres(cluster: str, amps):
    import pandas as pd
    import pickle as pkl
    from glob import glob

    
    df = pd.read_table(f'analysis/{cluster}_ANI.tsv',
               sep='\t',
               header=None,
               names=['s1', 's2', 'ANI', 'ortholog_frags', 'total_frags'])

    df['s1'] = df.s1.apply(lambda x: x.replace('.fna.gz', ''))
    df['s2'] = df.s2.apply(lambda x: x.replace('.fna.gz', ''))
    df['OF'] = df.ortholog_frags * 100 / df.total_frags
    df = df.drop(['ortholog_frags', 'total_frags'], axis=1)
    df['cluster'] = cluster
        
    b = pkl.load(open(f'analysis/{cluster}_seqs.pkl', 'rb'))
    print(f'loaded {len(b)} proteins')
    
    levels = pd.read_table(f'analysis/{cluster}.levels_assessment.tsv.gz',
                   sep='\t',
                   header='infer')
    levels = levels[['seq', 'L-III']]
    levels = levels.set_index('seq').to_dict()
    levels = levels['L-III']
    print(f'loaded {len(levels)} protein families')
    
    nlist = list()
    for _, s1, s2, og, of, _ in df.itertuples():
        if s1 != s2:
            s1genes, s1fams = (set(), set())
            s2genes, s2fams = (set(), set())
            for k, v in b.items():
                v = {x.split('_')[0] for x in v}
                family = levels.get(k, 'NA')
                if s1 in v:
                    s1genes.add(k)
                    s1fams.add(family)
                if s2 in v:
                    s2genes.add(k)
                    s2fams.add(family)
                
            total, shared = len(s1genes.union(s2genes)), len(s1genes.intersection(s2genes))
            ftotal, fshared = len(s1fams.union(s2fams)), len(s1fams.intersection(s2fams))
            
            a1 = set(amps.loc[amps['sample'] == s1, 'amp'])
            a2 = set(amps.loc[amps['sample'] == s2, 'amp'])
            tamp, sharedamp = len(a1.union(a2)), len(a1.intersection(a2))
            
            fa1 = set(amps.loc[amps['sample'] == s1, 'fam'])
            fa2 = set(amps.loc[amps['sample'] == s2, 'fam'])
            ftamp, fsharedamp = len(fa1.union(fa2)), len(fa1.intersection(fa2))
            
            nlist.append([s1, s2, cluster, og, of,
                          len(s1genes), len(s2genes),
                          total, shared,
                          len(s1fams), len(s2fams),
                          ftotal, fshared,
                          len(a1), len(a2),
                          tamp, sharedamp,
                          len(fa1), len(fa2),
                          ftamp, fsharedamp])
    
    newdf = pd.DataFrame(nlist,
                         columns=['sample1',
                                  'sample2',
                                  'cluster',
                                  'ANI',
                                  'OF',
                                  's1_proteins',
                                  's2_proteins',
                                  'total_nr_proteins',
                                  'shared_proteins',
                                  's1_protein_families',
                                  's2_protein_families',
                                  'total_protein_families',
                                  'shared_protein_families',
                                  's1_amps',
                                  's2_amps',
                                  'total_nr_amps',
                                  'shared_amps',
                                  's1_amp_families',
                                  's2_amp_families',
                                  'total_nr_amp_families',
                                  'shared_amp_families'])

    print(newdf.head())

    df.to_csv(f'analysis/{cluster}_ANI.tsv',
              sep='\t',
              header=None,
              index=None)

    newdf.to_csv(f'analysis/{cluster}_ANI_vs_proteinpep_overlap.tsv.gz', 
         sep='\t',
         header=True, 
         index=None)
         
