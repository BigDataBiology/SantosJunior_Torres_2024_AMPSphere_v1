def load_data():
    import pandas as pd
    lv1 = pd.read_table('analysis/output_clustering_significance_levelI.tsv')
    lv2 = pd.read_table('analysis/output_clustering_significance_levelII.tsv')
    lv3 = pd.read_table('analysis/output_clustering_significance_levelIII.tsv')
    return [lv1, lv2, lv3]
    

def plot_fig():
    import seaborn as sns
    import matplotlib.pyplot as plt
    from math import log
    
    cl = {1: 'I', 2: 'II', 3: 'III'}
    
    lv1, lv2, lv3 = load_data()
    lv1['min_id'] = lv1[['identity', 'gap_identity']].min(axis=1)
    lv2['min_id'] = lv2[['identity', 'gap_identity']].min(axis=1)
    lv3['min_id'] = lv3[['identity', 'gap_identity']].min(axis=1)

    lv1['Log(E-value)'] = lv1.evalue.apply(lambda x: log(x, 10))
    lv2['Log(E-value)'] = lv2.evalue.apply(lambda x: log(x, 10))
    lv3['Log(E-value)'] = lv3.evalue.apply(lambda x: log(x, 10))
    
    fig, axarr = plt.subplot_mosaic([['a)', 'b)', 'c)']], constrained_layout=True)
    sns.scatterplot(ax=axarr['a)'], data=lv1, x='min_id', y='Log(E-value)', hue='replicate', s=2.5, alpha=0.5)
    sns.scatterplot(ax=axarr['b)'], data=lv2, x='min_id', y='Log(E-value)', hue='replicate', s=2.5, alpha=0.5, legend=False)
    sns.scatterplot(ax=axarr['c)'], data=lv3, x='min_id', y='Log(E-value)', hue='replicate', s=2.5, alpha=0.5, legend=False)

    for idx, (label, ax) in enumerate(axarr.items()):
        ax.set_title(f'Clstr. Lv. {cl[idx+1]}', fontfamily='Sans Serif', fontstyle='italic')
        ax.axhline(y=log(1e-5, 10), color='black', linestyle='dashed', linewidth=1.0)
        ax.set_xlabel('Identity (%)')
        if idx > 0:
            ax.set_ylabel(None)

    plt.savefig('analysis/clustering_significance.svg')
    
plot_fig()

