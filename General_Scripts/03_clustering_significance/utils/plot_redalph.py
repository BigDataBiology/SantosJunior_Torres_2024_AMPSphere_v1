def load_data():
    import pandas as pd
    lI = pd.read_table('output_clustering_redalph_significance_levelI.tsv')
    lII = pd.read_table('output_clustering_redalph_significance_levelII.tsv')
    lIII = pd.read_table('output_clustering_redalph_significance_levelIII.tsv')
    return lI, lII, lIII
    
    
def process(dfin):
    import pandas as pd
    from collections import Counter
    from numpy import mean, std
    dfout = []
    for record in dfin.groupby('replicate'):
        s = Counter(record[1]['sig.'])
        s = s['*'] / sum(s.values())
        dfout.append(s)
    return mean(dfout), std(dfout)
    

def makedfs(lI, lII, lIII):
    import pandas as pd
    dfout = []
    dfout.append(process(lI))
    dfout.append(process(lII))
    dfout.append(process(lIII))
    df = pd.DataFrame(dfout,
                      columns=['mean', 'std'],
                      index=['level I', 'level II', 'level III'])
    return df
   

def siggraph(df):
    import matplotlib.pyplot as plt
    df['mean'].plot(kind='bar',
                    yerr=df['std'],
                    colormap='YlOrBr',
                    edgecolor='black',
                    grid=False)
    plt.ylabel('Proportion of significant alignments')
    plt.tight_layout()
    plt.savefig('redalph_significance_across_clustering_levels.svg')
    plt.close()
            
        
def idgraph(lI, lII, lIII):
    import matplotlib.pyplot as plt
    lI.identity.plot.density(label='level I', alpha=0.25)
    lII.identity.plot.density(label='level II', alpha=0.25)
    lIII.identity.plot.density(label='level III', alpha=0.25)
    plt.xlim(0, 100)
    plt.legend()
    plt.xlabel('Identity')
    plt.tight_layout()
    plt.savefig('redalph_identity_distribution.svg')
    plt.close()    

    
def covgraph(lI, lII, lIII):
    import matplotlib.pyplot as plt
    lI.coverage.plot.density(label='level I', alpha=0.25)
    lII.coverage.plot.density(label='level II', alpha=0.25)
    lIII.coverage.plot.density(label='level III', alpha=0.25)
    plt.xlim(0, 100)
    plt.legend()
    plt.xlabel('Coverage')
    plt.tight_layout()
    plt.savefig('redalph_coverage_distribution.svg')
    plt.close()

    
def redalph_graph():
    print('load info')
    lI, lII, lIII = load_data()
    print('convert into proportions')
    df = makedfs(lI, lII, lIII)
    print('plot significance')
    siggraph(df)
    print('plot identity')
    idgraph(lI, lII, lIII)
#    print('plot coverage')
#    covgraph(lI, lII, lIII)
    
