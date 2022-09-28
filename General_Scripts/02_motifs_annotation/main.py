def main():
    import sys
    import os
    
    from utils.annotation_complete import pipe_anno
    from utils.process_analysis import pa
    from utils.motif_plot import mplot
    from utils.comp_plot import compplot
    from utils.clustering_test import test_all
    
    os.makedirs('analysis', exist_ok=True)
    
    infile = 'data/AMPSphere_v.2022-03.faa.gz'
    ofile = 'AMPSphere_v.2022-03.annotation.tsv'
    
    if infile.endswith('gz'):
        import gzip
        infile = gzip.open(infile,
                           'rt',
                           encoding='utf-8')
    else:
        infile = open(infile, 'r')

    pipe_anno(infile, f'analysis/{ofile}')
    pa()
    mplot()
    test_all()
    compplot()


if __name__=='__main__':
    main()
    
