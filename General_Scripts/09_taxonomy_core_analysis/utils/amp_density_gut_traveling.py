def density_travel():
    import pandas as pd
    from scipy.stats import median_test
    
    metadata = "data/reduced_metadata.tsv"
    metadata = pd.read_table(metadata)
    metadata = metadata[metadata.general_envo_name.str.contains('human')]
    
    newdf = pd.DataFrame()
    for r in pd.read_table('data/bps-per-sample-per-taxon.tsv.xz',
                           header='infer',
                           sep='\t',
                           chunksize=5_000_000):
        r = r[r['sample'].isin(metadata.sample_accession)]
        r = r[r.level == 'species']
        r = r.drop(['taxid', 'level'], axis=1)
        newdf = pd.concat([newdf, r])
        print(len(newdf))
    
    newdf = newdf.merge(on='sample',
                        right=metadata[['sample_accession',
                                        'general_envo_name']].rename({'sample_accession': 'sample'}, axis=1))
    
    gmsc = pd.read_table('data/complete_amps_associated_taxonomy.tsv.gz')
    
    ampgenes = gmsc[gmsc.level == 'species'].groupby(['sample', 'source']).agg('size')
    ampgenes = ampgenes.reset_index().rename({'source': 'name', 0: 'amp_genes'}, axis=1)
    
    newdf = newdf.merge(on=['sample', 'name'], right=ampgenes)
    newdf['amp_density'] = newdf['amp_genes'] * 1e9 / newdf['nbps']
    
    habitats = ['human skin',
                'human mouth',
                'human respiratory tract',
                'human digestive tract',
                'human urogenital tract',
                'human gut']
    
    s = pd.DataFrame()
    for i in habitats:
        tdf = newdf[(newdf.nbps >= 1e6) & (newdf.general_envo_name == i)].sample(100)
        s = pd.concat([s, tdf])
    
    sns.boxplot(data=newdf[(newdf.nbps >= 1e6)],
                x='general_envo_name',
                y='amp_density',
                showfliers=False,
                color='white')
    
    sns.stripplot(data=s,
                  x='general_envo_name',
                  y='amp_density',
                  s=3)
                  
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig('analysis/species_amp_density_human_sites.svg')
    
    tdf = newdf[newdf.nbps >= 1e6].groupby("general_envo_name")['amp_density'].apply(lambda x: list(x))
    tdf = tdf.loc[habitats]
    
    stat, p, med, tbl = median_test(tdf.values[0],
                                    tdf.values[1],
                                    tdf.values[2],
                                    tdf.values[3],
                                    tdf.values[4],
                                    tdf.values[5])
    
    stat, p, med
    # (128.78401929137883, 4.311161122869419e-26, 761.6308649384755)
    
    k = newdf[newdf.nbps >= 1e6].groupby(['name', 'general_envo_name']).agg('size').sort_values().tail(10).index
    tdf = newdf[(newdf.nbps>=1e6) & (newdf.name.isin([x[0] for x in k])) & (newdf.general_envo_name == 'human gut')]
    
