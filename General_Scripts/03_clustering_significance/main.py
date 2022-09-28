from utils.clustering_analysis import cluster_analysis
from utils.plot_output import cluster_graph
from utils.clustering_redalph import cluster_redalph
from utils.plot_redalph import redalph_graph


def main():
    import os
    from os.path import exists
    
    os.makedirs('analysis/', exist_ok=True)
    
    cluster_analysis()
    cluster_graph()
    cluster_redalph()
    redalph_graph()
        
    ofiles = ['output_clustering_significance_levelIII.tsv',
              'output_clustering_significance_levelII.tsv',
              'output_clustering_significance_levelI.tsv',
              'output_clustering_redalph_significance_levelI.tsv',
              'output_clustering_redalph_significance_levelII.tsv',
              'output_clustering_redalph_significance_levelIII.tsv',
              'redalph_significance_across_clustering_levels.svg',
              'redalph_identity_distribution.svg',
              'redalph_coverage_distribution.svg',
              'coverage_distribution.svg',
              'significance_across_clustering_levels.svg',
              'identity_distribution.svg']
    
    for i in ofiles:
        if exists(i):
            os.rename(i,
                      f'analysis/{i}')
    
