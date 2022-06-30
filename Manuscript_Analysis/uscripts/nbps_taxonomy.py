def taxanalysis():
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    from itertools import combinations
    from scipy.stats import mannwhitneyu
    from statsmodels.stats.multitest import multipletests
    
    data = pd.read_table('data/amps_per_assembled_Gbp_per_taxon.tsv.xz')
    data.rename({'fixed': 'genus'}, axis=1, inplace=True)
    
    tax = pd.read_table("data/bac120_GTDB.tsv.xz",
                        sep='\t',
                        names=['domain', 'phylum',
                               'class', 'order',
                               'family', 'genus',
                               'species'])
                               
    tax.drop('species', axis=1, inplace=True)
    tax = tax.drop_duplicates()
    tax['genus'] = [x[3:] for x in tax.genus]
    
    data = data.merge(on='genus', right=tax)
    
    k = data.groupby('phylum').agg('size').sort_values().tail(13).index  # only phyla with 100 genera or more
    fdata = data[data.phylum.isin(k)].sort_values(by='phylum')
#    fdata = fdata[fdata.VAR_pct <= 100]
    porder = ['p__Planctomycetota', 'p__Bacteroidota', 'p__Chloroflexota',
              'p__Desulfobacterota', 'p__Verrucomicrobiota', 'p__Cyanobacteria',
              'p__Myxococcota', 'p__Acidobacteriota', 'p__Patescibacteria',
              'p__Proteobacteria', 'p__Actinobacteriota', 'p__Firmicutes',
              'p__Firmicutes_A']
    
    fig, ax = plt.subplots()
    sns.boxplot(ax=ax,
                data=fdata,
                x='phylum',
                y='amps_per_Gbp',
                showfliers=False,
                color='white',
                order=porder)
    
    sns.stripplot(ax=ax,
                  data=fdata,
                  x='phylum',
                  y='amps_per_Gbp',
                  s=2,
                  order=porder)
    
    ax.set_ylim(0, 2000)
    ax.set_xticklabels(porder, rotation=45)
    ax.set_xlabel('')
    ax.set_ylabel('AMPs per assembled gigabase pair')
    plt.tight_layout()
    fig.savefig('taxonomy.png', dpi=300)
    
    test = []
    for i, j in combinations(porder, 2):
        u, p = mannwhitneyu(fdata[fdata.phylum == i]['amps_per_Gbp'],
                            fdata[fdata.phylum == j]['amps_per_Gbp'])
        test.append((i, j, u, p))
    
    test = pd.DataFrame(test,
                        columns=['tax1', 'tax2',
                                 'u_stat', 'p-value'])
    
    _, test['p-value'], _, _ = multipletests(test['p-value'],
                                             method='bonferroni',
                                             is_sorted=False,
                                             returnsorted=False)
    
#    test = test[test['p-value'] < 5e-2]
    test.to_csv('taxonomy_test_ampspergbp.tsv',
                sep='\t',
                header=True,
                index=None)
                
