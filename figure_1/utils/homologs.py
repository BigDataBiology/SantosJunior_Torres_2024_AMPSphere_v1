def mmseqs(db, finput, output):
    '''
    Performs mmseqs searching
    '''
    import subprocess
    format=['query', 'target', 'evalue',
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
                    ','.join(format)])

                    
def search_homologs():
    '''
    Queries AMPSphere against different databases
    '''
    import os
    os.makedirs('homologs', exist_ok=True)
    outputs = {'DRAMP.fa': 'result_dramp.m8',
               'all_SmProt.fa.gz': 'result_SmProt.m8',
               'starPepDB.fasta': 'result_starPepDB.m8',
               'STsORFs.faa': 'result_STsORFs.m8'}
    for db, output in outputs.items():
        mmseqs(f'data/databases_homology/{db}',
               'data/AMPSphere_v.2021-03.faa.gz',
               f'homologs/{output}')


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


def search_GMGC():
    '''
    Performs the searching of AMPSphere against GMGC
    or in case user preferes, it can just skip this 
    process and retrieve the pre-computed resource
    '''
    import os
    import glob
    import lzma
    import subprocess
    import pandas as pd
    
    ask = input('''
    GMGC search takes much time and consumes much memory. 
    Because of that, we made available pre-computed results, if you prefer.
    In case you want to skip the homologs search in GMGC answer y,
    in the opposite case, reply with n.
    ''')
    
    format=['query', 'target', 'evalue',
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

        outdir = 'homologs/gmgc/'
        os.makedirs(outdir, exist_ok=True)

        for chunk in glob.glob('data/databases_homology/gmgc_chunks/*.fasta'):            
            outfile = chunk.split('/')[-1].replace('.fasta', '.m8')
            subprocess.call(['mmseqs',
                             'easy-search',
                             'data/AMPSphere_v.2021-03.faa.gz',
                             chunk,    
                             f'{outdir}/{outfile}',
                             'tmp',
                             '--disk-space-limit 5G',
                             '--db-load-mode 3',
                             '--threads 3',
                             '--format-output',
                             ','.join(format)])

        df = pd.DataFrame()
        for chunk in glob.glob('outdir/*.m8'):            
            dint = pd.DataFrame(chunk, sep='\t', header=None)
            dint = dint[dint[2] <= 1e-5]            
            df = pd.concat([df, dint])
            
        df.to_csv('homologs/truepep_gmgc_progenomes.m8',
                  sep='\t',
                  header=None,
                  index=None) 

        os.remove('homologs/gmgc')
        os.remove('data/databases_homology/gmgc_chunks')
        
    if (ask == 'Y') or (ask == 'y'):
        os.rename('data/databases_homology/truepep_gmgc_progenomes.m8.xz',
                  'homologs/gmgc.m8.xz')

        with lzma.open('homologs/gmgc.m8.xz', 'rt', encoding='utf-8') as infile:
            with open('homologs/gmgc.m8', 'w') as ofile:
                for row in infile:
                    ofile.write(row)

        os.remove('homologs/gmgc.m8.xz')


def compare_hits():
    '''
    Compare the hits obtained querying AMPSphere against different databases
    '''
    import glob
    import pandas as pd
    
    ofile = open('panelB_homologs_search.txt', 'a')

    all_candidates = set()
    for chunk in glob.glob('homologs/*.m8'):             
        homologs = pd.read_table(chunk, header=None)
        homologs = homologs[homologs[2] <= 1e-5][0]
        homologs = set(homologs)
        all_candidates = all_candidates.union(homologs)
        oname = chunk.split('/')[-1]
        oname = oname.replace('.m8', '')
        oname = oname.replace('result_', '')
        homologs = pd.DataFrame.from_dict(homologs)
        homologs.to_csv(f'homologs/{oname}_candidates.txt',
                        header=None,
                        index=None)
        print(f'It was detected in {oname}: {len(homologs)} AMPs with homologs',
              file=ofile)

    ofile.close()

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
    search_homologs()
    search_GMGC()
    all_candidates = compare_hits()
    overlaps(all_candidates)
        
