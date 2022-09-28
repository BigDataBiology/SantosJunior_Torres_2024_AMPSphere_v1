def test_clustering(infile, outfile, threshold):
    import pandas as pd
    from sklearn.cluster import AgglomerativeClustering  
    metrics = ['euclidean', 'braycurtis', 'correlation']
    methods = ['ward', 'complete', 'average', 'single']  
    data = pd.read_table(infile)
    if 'index' in data.columns: data = data.set_index('index')
    if 'Unnamed: 0' in data.columns: data = data.set_index('Unnamed: 0')
    clusters, techniques = [], []
    for me in metrics:
        for mo in methods:
            if mo == 'ward' and me != 'euclidean':
                pass
            else:
                for n in range(1, len(data)):
                    cmodel = AgglomerativeClustering(n_clusters=n,
                                                     affinity=me,
                                                     linkage=mo)
                    clusters.append(cmodel.fit_predict(data))
                    techniques.append([me, mo, n])
    techniques = ['_'.join([str(xi) for xi in x]) for x in techniques]
    df = pd.DataFrame(clusters,
                      columns=data.index,
                      index=techniques)
    df.to_csv(outfile, sep='\t', header=True, index=True)
    coords = []
    for idx_a, a in enumerate(df.columns):
        for idx_b, b in enumerate(df.columns):
            if idx_a < idx_b:
                 eq = len(df[df[a] == df[b]])
                 eq = (eq / len(df))
                 coords.append([a, b, eq])
    d = pd.DataFrame(coords, columns=['env1', 'env2', 'support'])
    print(d[d['support'] > threshold])


def test_all():
    print('Test top species')
    test_clustering('top_sps_clustermap.tsv', 'top_sps_cluster_analysis.tsv', 0.5)
    print('Test human body sites')
    test_clustering('human_sites_clustermap.tsv', 'human_sites_cluster_analysis.tsv', 0.5)
    print('Test guts')
    test_clustering('guts_clustermap.tsv', 'guts_cluster_analysis.tsv', 0.5)
    print('Test mega envo')
    test_clustering('megaenvo_clustermap.tsv', 'megaenvo_cluster_analysis.tsv', 0.5)

