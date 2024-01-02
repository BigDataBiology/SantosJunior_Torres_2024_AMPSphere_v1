from utils.singletons_handle import run_pipe
from utils.progenomes_amps import ampsphere2progenomes
from utils.features import calc_features
from utils.annotation import dramp_anno
from utils.cluster_analysis import process_cluster
from utils.metaG import metag
from utils.gmsc_genes import genes


def main():
    print('From prediction to families')
    run_pipe()
    print('Add progenomes genes')
    ampsphere2progenomes()
    print('Annotate metagenomes')
    metag()
    print('Generate fasta of genes')
    genes()
    print('Annotation with DRAMP')
    dramp_anno()
    print('Calculating AMP features')
    calc_features()
    print('Processing clusters')
    process_cluster()

if __name__ == '__main__':
    main()

