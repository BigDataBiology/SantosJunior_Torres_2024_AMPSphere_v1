def mfolder():
    import os
    os.makedirs('analysis/large_clusters', exist_ok=True)
    os.makedirs('analysis/large_clusters/aln', exist_ok=True)
    os.makedirs('analysis/large_clusters/rnacode', exist_ok=True)
    
    
def rnacode(in_file, out_file):
    import subprocess
    subprocess.call(['RNAcode',
                     '-o',
                     out_file,
                     '-t',
                     '--best-only',
                     '--stop-early',
                     '-n', '1000',
                     '-p', '0.05',
                     in_file])


def frag_rnacode():
    import os
    from glob import glob
    mfolder()
    for in_file in glob('*.aln'):
        ofile = in_file.replace('.aln', '.rnacode')
        rnacode(f'{in_file}',
                f'analysis/large_clusters/rnacode/{ofile}')
        os.rename(f'{in_file}',
                  f'analysis/large_clusters/aln/{in_file}')
        print(ofile)

