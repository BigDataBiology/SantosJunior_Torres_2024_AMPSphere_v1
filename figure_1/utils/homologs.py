def mmseqs(db, finput, output):
    '''
    Performs mmseqs searching
    '''
    import subprocess
    formlist=['query', 'target', 'evalue',
            'gapopen', 'pident', 'nident',
            'qstart', 'qend', 'qlen',
            'tstart', 'tend', 'tlen',
            'alnlen', 'raw', 'bits',
            'cigar', 'qseq', 'tseq',
            'qheader', 'theader', 'qaln',
            'taln', 'qframe', 'tframe',
            'mismatch', 'qcov', 'tcov']
    subprocess.call(['mmseqs',
                    'easy-search',
                    finput,
                    db,
                    output,
                    'tmp',
                    '--format-output',
                    ','.join(formlist)])

                    
def search_homologs():
    '''
    Queries AMPSphere against different databases
    '''
    import os
    os.makedirs('analysis/homologs', exist_ok=True)
    outputs = {'DRAMP.fa': 'result_dramp.m8',
               'all_SmProt.fa.gz': 'result_SmProt.m8',
               'starPepDB.fasta': 'result_starPepDB.m8',
               'STsORFs.faa': 'result_STsORFs.m8'}
    for db, output in outputs.items():
        mmseqs(f'data/databases_homology/{db}',
               'data/AMPSphere_v.2022-03.faa.gz',
               f'analysis/homologs/{output}')


def batch_iterator(iterator, batch_size):
    '''
    Returns lists of length batch_size.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.

    Available from: https://biopython.org/wiki/Split_large_file
    '''
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = iterator.next()
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch


def split_GMGC10():
    '''
    Split large GMGC database into smaller chunks
    '''
    from Bio import SeqIO
    import os

    gmgc_address = 'data/databases_homology/GMGC10.proGenomes.faa'
    outdir = 'data/databases_homology/gmgc_chunks/'

    os.makedirs(outdir, exist_ok=True)       

    record_iter = SeqIO.parse(open(gmgc_address), 'fasta')
    for i, batch in enumerate(batch_iterator(record_iter, 1_000_000)):
        filename = f'group_{i+1}.fasta'
        with open(f'{outdir}/{filename}', 'w') as handle:
            count = SeqIO.write(batch, handle, 'fasta')
        print(f'Wrote {count} records to {filename}')


def search_GMGC(ask):
    '''
    Performs the searching of AMPSphere against GMGC
    or in case user preferes, it can just skip this 
    process and retrieve the pre-computed resource

    :input: ask - str. Y or N to use the pre-computed
                  resource
    '''
    import os
    import glob
    import lzma
    import subprocess
    import pandas as pd
       
    formlist=['query', 'target', 'evalue',
        'gapopen', 'pident', 'nident',
        'qstart', 'qend', 'qlen',
        'tstart', 'tend', 'tlen',
        'alnlen', 'raw', 'bits',
        'cigar', 'qseq', 'tseq',
        'qheader', 'theader', 'qaln',
        'taln', 'qframe', 'tframe',
        'mismatch', 'qcov', 'tcov']

    if (ask == 'N') or (ask == 'n'):
        print('Splitting the large GMGC file into small chunks')
        split_GMGC10()

        outdir = 'analysis/homologs/gmgc/'
        os.makedirs(outdir, exist_ok=True)

        for chunk in glob.glob('data/databases_homology/gmgc_chunks/*.fasta'):            
            outfile = chunk.split('/')[-1].replace('.fasta', '.m8')
            subprocess.call(['mmseqs',
                             'easy-search',
                             'data/AMPSphere_v.2022-03.faa.gz',
                             chunk,    
                             f'{outdir}/{outfile}',
                             'tmp',
                             '--disk-space-limit 5G',
                             '--db-load-mode 3',
                             '--threads 3',
                             '--format-output',
                             ','.join(formlist)])

        df = pd.DataFrame()
        for chunk in glob.glob(f'{outdir}/*.m8'):            
            dint = pd.DataFrame(chunk, sep='\t', header=None)
            dint = dint[dint[2] <= 1e-5]            
            df = pd.concat([df, dint])
            
        df.to_csv('analysis/homologs/result_gmgc.m8',
                  sep='\t',
                  header=None,
                  index=None) 

        os.remove('analysis/homologs/gmgc')
        os.remove('data/databases_homology/gmgc_chunks')
        
    if (ask == 'Y') or (ask == 'y'):
        with  lzma.open('data/databases_homology/true_pep_2022_vs_progenomesgmgc.tsv.xz',
                        'rt',
                        encoding='utf-8') as infile:
            with open('analysis/homologs/result_gmgc.m8', 'w') as ofile:
                for row in infile:
                    ofile.write(row)


def compare_hits():
    '''
    Compare the hits obtained querying AMPSphere against different databases
    '''
    import glob
    import pandas as pd
    
    ofile = open('panelB_homologs_search.txt', 'a')

    all_candidates = set()
    for chunk in glob.glob('analysis/homologs/*.m8'):             
        homologs = pd.read_table(chunk, header=None)
        homologs = homologs[homologs[2] <= 1e-5][0]
        homologs = set(homologs)
        all_candidates = all_candidates.union(homologs)
        oname = chunk.split('/')[-1]
        oname = oname.replace('.m8', '')
        oname = oname.replace('result_', '')
        homologs = pd.DataFrame.from_dict(homologs)
        homologs.to_csv(f'analysis/homologs/{oname}_candidates.txt',
                        header=None,
                        index=None)
        print(f'It was detected in {oname}: {len(homologs)} AMPs with homologs',
              file=ofile)

    ofile.close()

    homologs = pd.DataFrame.from_dict(all_candidates)
    homologs.to_csv(f'analysis/homologs/annotated_candidates.txt',
                    header=None,
                    index=None)

    return all_candidates


def overlaps(all_candidates):
    '''
    Detects the overlap between quality candidates and all AMPs with homologs
    '''
    import pandas as pd
  
    ofile = open('panelB_homologs_search.txt', 'a')
    
    quality_candidates = set()
    infiles = ['high_quality_candidates.txt', 'quality_candidates.txt']
    for i in infiles:
        df = pd.read_table(f'data/{i}', header=None)
        df = set(df[0])
        quality_candidates = quality_candidates.union(df)
    
    intersect = all_candidates.intersection(quality_candidates)
    print(f'AMP candidates of quality with homologs: {len(intersect)}',
          file=ofile)

    QA = quality_candidates - all_candidates
    print(f'AMP candidates of quality without homologs: {len(QA)}',
          file=ofile)

    AQ = all_candidates - quality_candidates
    print(f'AMP with homologs failing quality tests: {len(AQ)}',
          file=ofile)

    ofile.close()
    
    
def homologs():
    import os
    from .timeout_input import timeout_input
    
    
    print('Diverse db searching')
    search_homologs()
    print('')
    t, answ = timeout_input('''
    GMGC search takes much time and consumes much memory. 
    Because of that, we made available pre-computed results, if you prefer.
    In case you want to skip the homologs search in GMGC answer y,
    in the opposite case, reply with n.''', 5, 'y')
    search_GMGC(answ)
    print('Getting all candidates')
    all_candidates = compare_hits()
    print('Generating numbers for Venns Diagram')
    overlaps(all_candidates)
    os.remove('tmp/')
        
