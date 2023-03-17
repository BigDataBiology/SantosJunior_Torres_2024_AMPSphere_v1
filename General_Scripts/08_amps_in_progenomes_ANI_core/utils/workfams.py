import os
import pickle
import itertools
import subprocess
import pandas as pd
import pickle as pkl
from glob import glob
from os.path import exists
from itertools import chain


def redseq(seq:str) -> str:
    '''
    # Clustering after reducing alphabet
    to 8 letters:
    [LVMIC], [AG], [ST], P,
    [FWY], [EDQN], [KR], H
    # Refs:
    # doi: 10.1093/protein/13.3.149
    # doi: 10.1093/bioinformatics/btp164
    # Solis, AD. Proteins 2015; 83:2198–2216.
    '''
    aas = {'V': 'L', 'M': 'L',
           'I': 'L', 'C': 'L',
           'G': 'A', 'T': 'S',
           'W': 'F', 'Y': 'F',
           'D': 'E', 'Q': 'E',
           'N': 'E', 'R': 'K'}
    nseq = ''
    for i in seq.strip():
        nseq += aas.get(i, i)       
    return nseq


def load_dict(f):
    return pickle.load(open(f, 'rb'))


def is_command(cmds):
    """Given a command returns its path, or None.
    Given a list of commands returns the first recoverable path, or None.
    """
    try:
        from shutil import which as which  # python3 only
    except ImportError:
        from distutils.spawn import find_executable as which
    if isinstance(cmds, str):
        return which(cmds)
    else:
        for cmd in cmds:
            path = which(cmd)
            if path is not None:
                return path
        return path


def call_cdhit(fin, fout, threshold, wordsize, threads, stdout=None):
    """
    Function to perform peptides clustering using cd-hit. Fixed
    parameters include: global identity, print alignment overlap,
    cluster sorting by decreasing size in the clstr and fasta
    files, slow clustering mode, band_width (5), throw away
    sequences under 6 -- safe margin, minimum coverage of 80% of
    the shorter sequence, print entire access code until first
    space, maximum disk memory usage of 16GB.
    """
    cdhit_exe = is_command(['cd-hit', 'cdhit'])
    if cdhit_exe is None:
        print('[ERROR] -- Cdhit Not Found --')
    subprocess.check_call(
        [cdhit_exe,
         '-i', fin,
         '-o', fout,
         '-c', str(threshold),
         '-n', str(wordsize),
         '-G', '1',
         '-g', '1',
         '-b', '5',
         '-l', '6',
         '-p', '1',
         '-sf', '1',
         '-sc', '1',
         '-aS', '0.8',
         '-M', '16000',
         '-d', '0',
         '-T', str(threads)],
        stdout=stdout)


def call_crev(fin, fout, stdout=None):
    """
    Function to combine a .clstr file with its parent .clstr file
    after a cd-hit hierarchical clustering procedure. It uses the
    script clstr_rev.pl in the PATH as clstr_rev.
    """
    crev_exe = is_command(['clstr_rev'])
    if crev_exe is None:
        print('[ERROR] -- Clstr_rev Not Found --')
    subprocess.check_call([crev_exe, fin, fout], stdout=stdout)


def delimit(iterable, splitstring):
    """
    Function to split cd-hit cluster files in a sensible way.
    Groups the items in the list using a keyword, and creates a
    list of sublists that are grouped accordingly.
    Returns: a list of sublists.
    """
    return [list(g) for k, g in itertools.groupby(
        iterable, lambda x:x in splitstring) if not k]


def parse_clusters(l, cluster_file, reference, foutput):
    """
    Reads a cluster file and generates a list of the contents, where
    the string Cluster denotes the beginning of a new sequence cluster.
    The delimit() function is used to create a list of sublists based on the
    appearance of the string Cluster in the original list. The list of sublists
    is reverse sorted and is written as a new sorted clstr file.
    l = level (I,II,III)
    level depends on cluster file:
    level I -- nr100
         II -- nr100-85
        III -- nr100-85-75
    foutput = clusters sorted by size with AMPL headers reference
    """
    resource = pd.read_table(reference, sep='\t', header='infer')
    resource = resource[['sequence', 'accession']]
    resource = resource.set_index('accession')
    resource = resource.to_dict()['sequence']
    lines = []
    with open(cluster_file, 'r') as fh:
        for line in fh:
            if line.startswith('>'):
                lines.append("Cluster")
            else:
                lines.append(line.strip())
    sublists = delimit(lines, ("Cluster",))
    sublists.sort(key=len, reverse=True)
    k = ['... at ', '... ', ', >', ' ']
    with open(foutput, 'w') as db:
        for i in range(len(sublists)):
            for j in range(len(sublists[i])):
                cl_row = sublists[i][j]
                for kc in k: cl_row = cl_row.replace(kc, '\t')
                cl_row = cl_row.split('\t')
                hseq = cl_row[2]
                fmil = str(i).zfill(6)
                fmil = '_'.join(fmil[m:m + 3] for m in range(0, len(fmil), 3))
                db.writelines(f'L-{l}.{fmil}\t{resource[hseq]}\t{cl_row[2]}\t{cl_row[1]}\t{cl_row[3]}\n')


def editdictwoclones(cluster):
    '''
    Keeps only one clone per cluster inside
    a species, reducing the dictionary and
    the data redundancy
    '''   
    print(f'WORKING on cluster {cluster}')
    record = f'ani_amp_res/{cluster}_strain_clones.tsv'
    
    if exists(record):
        df = pd.read_table(record, index_col=0)
    else: 
        return False
        
    # load proteins
    x = load_dict(f'analysis/{cluster}_seqs.pkl')
    samples = set()
    for w in x.values():
        for wx in w:
            samples.add(wx.split('_')[0])
    
    # generating list of genomes to del from first
    samples_ani = set(df.index)
    todel2 = samples_ani - samples
    if len(todel2) > 0:
        df = df.drop(todel2, axis=0) 
        df = df.groupby('clone').apply(lambda x: x.reset_index().sample(1))
        df = df.reset_index(drop=True)['index'].tolist()
    else:
        df = list(df.index)

    # prunning protdict
    if len(df) > 10:
        ofile = open('data_from_clusters_wo_clones.txt', 'a')
        print(f'After eliminating clones\n\tthere is a total of {len(df)} genomes for this cluster')
        x = {k: [w for w in v if w.split('_')[0] in df] for k, v in x.items()}
        x = {k: v for k, v in x.items() if len(v) > 0}
        print(f'\t{len(x)} proteins were left in this cluster')
        ofile.write(f'{cluster}\t{len(df)}\t{len(x)}\n')
        ofile.close()
        return x
        
    else:
        return False
    
    
def cluster_fams(cluster, Nthreads=3, verbose=False):
    '''
    Use the same criteria used to cluster AMP families
    and the same method based in the CDHit, we create
    then a table with all protein clusters, later used
    to get the number of core, shell, and accessory
    protein families per each prevalence cutoff per
    species
    '''
#    if exists(f'analysis/{cluster}_families.tsv.gz'):
#        return f'... {cluster}'

    b = editdictwoclones(cluster)

    if not b:
        return (False, False)

    with open('intfasta.tmp.fa', 'wt') as ofile:
        for idx, (k, v) in enumerate(b.items()):
            ofile.write(f'>T_{idx}\n{redseq(k)}\n')

    if verbose: print('\t\t1st stage - 100% of identity')
    with open('nr100.clstr.log', 'w') as log:
        call_cdhit(
            fin='intfasta.tmp.fa',
            fout='nr100',
            threshold=1,
            wordsize=5,
            threads=Nthreads,
            stdout=log)

    if verbose: print('\t\t2nd stage -  85% of identity')
    with open(f'nr85.clstr.log', 'w') as log:
        call_cdhit(
            fin='nr100',
            fout='nr85',
            threshold=0.85,
            wordsize=5,
            threads=Nthreads,
            stdout=log)

    if verbose: print('\t\t3rd stage -  75% of identity')
    with open('nr75.clstr.log', 'w') as log:
        call_cdhit(
            fin='nr85',
            fout='nr75',
            threshold=0.75,
            wordsize=5,
            threads=Nthreads,
            stdout=log)

    if verbose: print('Converting output files to a human readable...')
    with open('nr100-85.clstr', 'w') as out:
        call_crev('nr100.clstr',
                  'nr85.clstr',
                  out)
                  
    with open('nr100-85-75.clstr', 'w') as out:
        call_crev('nr100-85.clstr',
                  'nr75.clstr',
                  out)
    
    hs = []
    for idx, (k, v) in enumerate(b.items()):
        hs.append((f'T_{idx}', k, v))
        
    hs = pd.DataFrame(hs,
                      columns=['accession', 'sequence', 'samples'])

    hs.to_csv('ref.tsv',
              sep='\t',
              header=True,
              index=None)

    if verbose: print('Parsing clusters from level I')
    parse_clusters('I',
                   f'nr100.clstr',
                   f'ref.tsv',
                   f'nr100.sorted.clstr')

    if verbose: print('Parsing clusters from level II')
    parse_clusters('II',
                   'nr100-85.clstr',
                   'ref.tsv',
                   'nr100-85.sorted.clstr')

    if verbose: print('Parsing clusters from level III')
    parse_clusters('III',
                   'nr100-85-75.clstr',
                   'ref.tsv',
                   'nr100-85-75.sorted.clstr')
    
    df_100 = pd.read_table(f'nr100.sorted.clstr',
                           sep='\t',
                           header=None,
                           names=['L-I', 'seq',
                                  'sequence', 'length',
                                  'clustering'])

    df_85 = pd.read_table(f'nr100-85.sorted.clstr',
                           sep='\t',
                           header=None,
                           names=['L-II', 'seq',
                                  'sequence', 'length',
                                  'clustering'])

    df_75 = pd.read_table(f'nr100-85-75.sorted.clstr',
                           sep='\t',
                           header=None,
                           names=['L-III', 'seq',
                                  'sequence', 'length',
                                  'clustering'])
    
    df_100 = df_100[['seq', 'sequence', 'L-I']]
    df_85 = df_85[['seq', 'L-II']]
    df_75 = df_75[['seq', 'L-III', 'clustering']]
    
    data = df_100.merge(on='seq', right=df_85)
    data = data.merge(on='seq', right=df_75)
    
    k = data['L-III'].value_counts()
    k = k[k>=8]
    
    data = data[data['L-III'].isin(k.index)]

    data.to_csv(f'analysis/{cluster}.levels_assessment.tsv.gz',
                sep='\t',
                header=True,
                index=None)

    return data, b
    
    
def workfam(f, verbose=False):
    cluster = f.split('/')[-1]
    cluster = cluster.split('_seqs')[0]
    data, b = cluster_fams(cluster, Nthreads=3, verbose=verbose)
    
    if isinstance(data, bool):
        return f'ERROR {cluster} not included in ANI'
    else:            
        data = data.loc[:, ['seq', 'L-III']]
        data['samples'] = data['seq'].apply(lambda x: b[x])   
        data = data.drop('seq',
                         axis=1).groupby('L-III').apply(lambda x: list(chain.from_iterable(x['samples'])))        
      
        data = data.apply(lambda x: set([j.split("_")[0] for j in x]))
        data = data.reset_index().rename({0: 'samples'}, axis=1)                 
        data['length'] = data.samples.apply(lambda x: len(x))
    
        N = len({x.split('_')[0] for x in list(chain.from_iterable(b.values()))})
        data['perc'] = data['length'] * 100 / N
    
        fams = []
        for c in range(0, 101, 5):
           n = len(data[data.perc >= c])
           p = n*100 / len(data)
           fams.append((c, n, p))
        
        fams = pd.DataFrame(fams,
                            columns=['cutoff_prevalence',
                                     'number of families',
                                     'percent of families'])

        fams.to_csv(f'analysis/{cluster}_families.tsv.gz',
                    sep='\t',
                    header=True,
                    index=None)
    
        if verbose: print('Organizing results')
        flist = ['nr75', 'nr85', 'nr100',
                 'nr75.clstr', 'nr85.clstr', 'nr100.clstr',
                 'nr100-85.clstr', 'nr100-85-75.clstr',
                 'nr100-85-75.sorted.clstr', 'nr100-85.sorted.clstr',
                 'nr100.sorted.clstr', 'nr75.clstr.log',
                 'nr85.clstr.log', 'nr100.clstr.log',
                 'intfasta.tmp.fa', 'ref.tsv']
             
        for x in flist: os.remove(f'{x}')         
    
    
def families():
    from glob import glob
    for infile in glob('analysis/*_seqs.pkl'): 
        workfam(infile,
                verbose=False)


if __name__ == '__main__':
    families()
    
