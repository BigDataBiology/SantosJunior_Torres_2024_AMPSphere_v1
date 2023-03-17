def eliminate_non_standard_aas(infile):
    '''
    Eliminate rows with non-standard residues and
    check for the consistency of the original
    prediction results table

    Input: 
    infile - address of file containing the output from macrel
             prediction for metagenomes
    
    Output:
    ofile - address of file to be saved containing filtered results
    '''
    import gzip
    import re
    
    amps_seen = []
    with gzip.open(infile, 'rt', encoding='utf-8') as f:
        for line in f:
            line = line.split('\t')
            if len(line) == 6: 
                matches = re.findall('b|B|o|O|x|X|u|U|z|Z|j|J', line[1])
                if matches: continue
                else: amps_seen.append(line[1])

    print(f'It was found a total of {len(amps_seen)} AMPs')
    
    return amps_seen
    

def out_singletons(amps_seen):
    '''
    Eliminating AMPs do not appearing more than once
    '''
    from collections import Counter

    counts = Counter(amps_seen)
    singletons, non_singletons = set(), dict()
    for k, v in counts.items():
        if v == 1: singletons.add(k) 
        else: non_singletons[k] = v

    print(f'Singletons: {len(singletons)}')
    print(f'Non-Singletons: {len(non_singletons)}')

    return (singletons, non_singletons)        
    
    
def precomputed_res(data_folder, outdir):
    '''
    Uploads and returns the singleton candidates
    saved by alignment against DRAMP using an old
    version of MMSeqs (same as in their original 
    publication)
    '''
    import pandas as pd
    
    data = pd.read_table(f'{data_folder}/DRAMP_filter.raw.tsv.xz',
                         sep='\t',
                         header=None,
                         names=['query', 'target', 'fident',
                                'alnlen', 'mismatch', 'gapopen',
                                'qstart', 'qend', 'tstart', 'tend',
                                'evalue', 'bits'])

    parsed = data[(data.fident >= 0.75) & (data.evalue <= 1e-5)]

    parsed = parsed.sort_values(by=['bits', 'evalue', 'fident'],
                                ascending=[False, True, False])
      
    parsed.to_csv(f'{outdir}/DRAMP_filter.parsed.tsv.xz', sep='\t', header=True, index=None)
    singletons_saved = set(parsed['query'])

    print(f'We could save {len(singletons_saved)}')

    return singletons_saved


def recovering_singletons(singletons, outdir, dramp):
    '''
    Recover singletons with homologs in DRAMP,
    a database of AMPs
    '''
    import lzma
    import tempfile
    import pandas as pd
    from .utils import call_mmseqs
     
    print ('Creating temporary fasta of singletons to search for DRAMP homologs')
    with tempfile.TemporaryDirectory() as tmpdirname:
        with open(f'{tmpdirname}/temp.fa', 'w') as tf:
            for idx, s in enumerate(singletons):
                tf.write(f'>{idx}\n{s}\n')

        print ('Searching reminiscent singletons against DRAMP, it can take a while')
        bout = f'{tmpdirname}/DRAMP_filter.raw.tsv'
        fin = f'{tmpdirname}/temp.fa'
        Nthreads = 3
        call_mmseqs(fin, bout, dramp, Nthreads, tmp='tmp', stdout=None)
        
        data = pd.read_table(bout,
                             sep='\t',
                             header=None,
                             names=['query', 'target', 'fident',
                                    'alnlen', 'mismatch', 'gapopen',
                                    'qstart', 'qend', 'tstart', 'tend',
                                    'evalue', 'bits'])
        
        data.to_csv(bout+'.xz', sep='\t', header=True, index=None)
        
        parsed = data[(data.fident >= 0.75) & (data.evalue <= 1e-5)]
        parsed = parsed.sort_values(by=['bits', 'evalue', 'fident'],
                                    ascending=[False, True, False])
               
        singletons = list(singletons)
        parsed['query'] = [singletons[x] for x in parsed['query']]

        parsed.to_csv(f'{outdir}/DRAMP_filter.parsed.tsv.xz', sep='\t', header=True, index=None)
        singletons_saved = set(parsed['query'])

        print(f'We could save {len(singletons_saved)}')

        return singletons_saved

    
def reduceseq(seq):
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


def clusteramps(non_singletons, outdir):
    '''
    Recover singletons with homologs in DRAMP,
    a database of AMPs
    '''
    import os
    from .utils import call_cdhit, call_crev

    Nthreads = 3
    
    with open(f'{outdir}/temp.fa', 'w') as infile:
        for k in non_singletons.keys():
            infile.write(f'>{k}\n{reduceseq(k)}\n')

    print('Hierarchical clustering started')

    print('\t\t1st stage - 100% of identity')
    with open(f'{outdir}/nr100.clstr.log', 'w') as log:
        call_cdhit(
            fin=f'{outdir}/temp.fa',
            fout=f'{outdir}/nr100',
            threshold=1,
            wordsize=5,
            threads=Nthreads,
            stdout=log)

    print('\t\t2nd stage -  85% of identity')
    with open(f'{outdir}/nr85.clstr.log', 'w') as log:
        call_cdhit(
            fin=f'{outdir}/nr100',
            fout=f'{outdir}/nr85',
            threshold=0.85,
            wordsize=5,
            threads=Nthreads,
            stdout=log)

    print('\t\t3rd stage -  75% of identity')
    with open(f'{outdir}/nr75.clstr.log', 'w') as log:
        call_cdhit(
            fin=f'{outdir}/nr85',
            fout=f'{outdir}/nr75',
            threshold=0.75,
            wordsize=5,
            threads=Nthreads,
            stdout=log)

    print('Converting output files to a human readable...')

    with open(f'{outdir}/nr100-85.clstr', 'w') as out:
        call_crev(f'{outdir}/nr100.clstr',
                  f'{outdir}/nr85.clstr',
                  out)

    with open(f'{outdir}/nr100-85-75.clstr', 'w') as out:
        call_crev(f'{outdir}/nr100-85.clstr',
                  f'{outdir}/nr75.clstr',
                  out)

    print('Cleaning environment')
    os.remove(f'{outdir}/temp.fa')


def AMPSPHERE_code(idx):
    '''
    Generate a code access for AMPSphere
    from an index number for the sorted
    input of AMPs by their redundancy size
    '''
    ac = str(idx).zfill(6)
    ac = '_'.join([ac[idx:idx + 3]
                   for idx, val in enumerate(ac) if idx % 3 == 0])

    return f'AMP10.{ac}'

    
def link_gmsc(infile2, datafile, outfile):
    '''
    Function to generate the correspondence list between the
    headers from AMPsphere and from GMSC
    '''
    import gzip
    import pandas as pd
    from collections import defaultdict
        
    # get descending sorted amps list
    amp_list = pd.DataFrame.from_dict(infile2, orient='index', columns=['size'])
    amp_list = amp_list.reset_index()
    amp_list = amp_list.rename({'index': 'sequence'}, axis=1)
    amp_list['length'] = amp_list['sequence'].str.len()

    amp_list = amp_list.sort_values(by=['size',
                                        'length',
                                        'sequence'],
                                    ascending=[False,
                                               True,
                                               True])
    
    amp_list = amp_list.reset_index(drop=True)                                           
    amp_list['accession'] = [AMPSPHERE_code(idx) for idx in amp_list.index]

    # defining references
    print('\t\tOpening initial files, it can take a while...')
    seqs = set(amp_list['sequence'])                
    dictresource = defaultdict(list)
    with gzip.open(datafile, 'rt', encoding='utf-8') as enc:
        for x in enc:
            tokens = x.split('\t')
            if tokens[1] in seqs:
                dictresource[tokens[1]].append(tokens[0])
    
    dictresource = [(k, ','.join(v), len(v)) for k,v in dictresource.items()]
    
    genes = pd.DataFrame(dictresource,
                         columns=['sequence',
                                  'genes',
                                  'n_of_genes'])
      
    df = amp_list.merge(on='sequence', right=genes)

    # checking if genes and samples are the same
    df = df[df['size'] == df['n_of_genes']]

    # reordering the df
    df = df[['accession', 'sequence', 'length', 'n_of_genes', 'genes']]
    df.to_csv(outfile, sep='\t', header=True, index=None)


def delimit(iterable, splitstring):
    """
    Function to split cd-hit cluster files in a sensible way.
    Groups the items in the list using a keyword, and creates a
    list of sublists that are grouped accordingly.
    Returns: a list of sublists.
    """
    import itertools
    
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

    foutput = clusters sorted by size with AMPsphere headers reference
    """
    import pandas as pd
    
    resource = pd.read_table(reference, sep='\t', header='infer')
    resource = resource[['sequence', 'accession']]
    resource = resource.set_index('sequence')
    resource = resource.to_dict()['accession']
    
    # Collect the file contents then delimit lists based on presence of
    # 'Cluster' in text
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
                db.writelines(f'SPHERE-{l}.{fmil}\t{resource[hseq]}\t{cl_row[2]}\t{cl_row[1]}\t{cl_row[3]}\n')


def spheres(analysis_folder):
    '''
    Transform the scattered info of clustering
    into an unified table
    '''
    import os
    import pandas as pd
    
    print('Parsing clusters from level I')
    parse_clusters('I',
                   f'{analysis_folder}/nr100.clstr',
                   f'{analysis_folder}/AMPsphere_GMSC_correspondence.tsv.gz',
                   f'{analysis_folder}/nr100.sorted.clstr')

    print('Parsing clusters from level II')
    parse_clusters('II',
                   f'{analysis_folder}/nr100-85.clstr',
                   f'{analysis_folder}/AMPsphere_GMSC_correspondence.tsv.gz',
                   f'{analysis_folder}/nr100-85.sorted.clstr')

    print('Parsing clusters from level III')
    parse_clusters('III',
                   f'{analysis_folder}/nr100-85-75.clstr',
                   f'{analysis_folder}/AMPsphere_GMSC_correspondence.tsv.gz',
                   f'{analysis_folder}/nr100-85-75.sorted.clstr')
    
    print('Parsing clusters together')

    df_100 = pd.read_table(f'{analysis_folder}/nr100.sorted.clstr',
                           sep='\t',
                           header=None,
                           names=['SPHERE-I', 'amp',
                                  'sequence', 'length',
                                  'clustering'])

    df_85 = pd.read_table(f'{analysis_folder}/nr100-85.sorted.clstr',
                           sep='\t',
                           header=None,
                           names=['SPHERE-II', 'amp',
                                  'sequence', 'length',
                                  'clustering'])

    df_75 = pd.read_table(f'{analysis_folder}/nr100-85-75.sorted.clstr',
                           sep='\t',
                           header=None,
                           names=['SPHERE-III', 'amp',
                                  'sequence', 'length',
                                  'clustering'])

    df_100 = df_100[['amp', 'sequence', 'SPHERE-I']]
    df_85 = df_85[['amp', 'SPHERE-II']]
    df_75 = df_75[['amp', 'SPHERE-III', 'clustering']]

    data = df_100.merge(on='amp', right=df_85)
    data = data.merge(on='amp', right=df_75)

    data = data.rename({'amp': 'AMP accession',
                        'clustering': 'evaluation vs. representative',
                        'SPHERE-I': 'SPHERE_fam level I',
                        'SPHERE-II': 'SPHERE_fam level II',
                        'SPHERE-III': 'SPHERE_fam level III'},
                       axis=1)

    data = data[['AMP accession',
                 'evaluation vs. representative',
                 'SPHERE_fam level I', 'SPHERE_fam level II',
                 'SPHERE_fam level III']]

    data = data.sort_values(by='AMP accession')             
    data.to_csv(f'{analysis_folder}/SPHERE_v.2022-03.levels_assessment.tsv.gz',
                sep='\t', header=True, index=None)
                
    print('Organizing results')

    flist = ['nr75', 'nr85', 'nr100',
             'nr75.clstr', 'nr85.clstr', 'nr100.clstr',
             'nr100-85.clstr', 'nr100-85-75.clstr',
             'nr100-85-75.sorted.clstr', 'nr100-85.sorted.clstr',
             'nr100.sorted.clstr', 'nr75.clstr.log',
             'nr85.clstr.log', 'nr100.clstr.log']

    os.makedirs(f'{analysis_folder}/clustering', exist_ok=True)

    for f in flist: os.replace(f'{analysis_folder}/{f}',
                               f'{analysis_folder}/clustering/{f}')


def fasta_file(analysis_folder):
    import gzip
    import pandas as pd
    
    data = pd.read_table(f'{analysis_folder}/AMPsphere_GMSC_correspondence.tsv.gz',
                         sep='\t', header='infer')
    
    fams = pd.read_table(f'{analysis_folder}/SPHERE_v.2022-03.levels_assessment.tsv.gz',
                         sep='\t', header='infer')
    
    fams = fams.rename({'AMP accession': 'accession',
                        'evaluation vs. representative': 'clustering',
                        'SPHERE_fam level I': 'SPHERE-I',
                        'SPHERE_fam level II': 'SPHERE-II',
                        'SPHERE_fam level III': 'SPHERE-III'},
                       axis=1)
    
    df = data[['accession', 'sequence']]
    df = df.merge(on='accession', right=fams[['accession', 'SPHERE-III', 'clustering']])
    
    with gzip.open(f'{analysis_folder}/AMPSphere_v.2022-03.faa.gz', 'wt', encoding='utf-8') as ofile:
        for row in df.itertuples():
            _, header, seq, f, _ = row
            ofile.write(f'>{header} | {f}\n{seq}\n')

    print('Generating AMPsphere only with the representatives...')

    with open(f'{analysis_folder}/AMPsphere_representatives.faa', 'w+') as ofile:
        for row in df.itertuples():
            _, header, seq, f, clustering = row
            if clustering == '*':
                ofile.write(f'>{header} | {f}\n{seq}\n')
   
               
def run_pipe():
    '''
    Run the entire pipeline until the clustering
    and organization of families
    '''
    import os
    from .timeout_input import timeout_input
    
    print('Set up environment')
    data_folder = 'data/'
    analysis_folder = 'analysis/'

    dramp = f'{data_folder}/DRAMP_GCP2019.fasta'

    for d in [data_folder, analysis_folder]:        
        os.makedirs(d, exist_ok=True)

    print('Eliminating sequences containing non-standard residues')
    amps_seen = eliminate_non_standard_aas(f'{data_folder}/GMSC10.Macrel_05.AMPs.tsv.gz')
    
    print('-- Filtering singletons:')
    singletons, non_singletons = out_singletons(amps_seen)
    
    print('Recovering AMPs matching to DRAMP but still singletons')
    
    t, answ = timeout_input('''There was a change in the versions of MMSeqs2.
                 This impacted the results in this step. Do you want to 
                 proceed the alignment? Alternatively, there is a pre-computed
                 set of results which can be used.
                 ''', 5, 'n')
                 
    if answ in ['y', 'Y', 'yes', 'Yes', 'YES']:
        singletons_saved = recovering_singletons(singletons,
                                                 analysis_folder,
                                                 dramp)
    else:
        print('Assuming the precomputed results')        
        singletons_saved = precomputed_res(data_folder,
                                           analysis_folder)

    for s in singletons_saved: non_singletons[s] = 1
     
    print('Clustering AMPs')
    clusteramps(non_singletons, analysis_folder)
    
    print('Link AMPs to GMSC genes')
    link_gmsc(non_singletons,
              f'{data_folder}/GMSC10.Macrel_05.AMPs.tsv.gz',
              f'{analysis_folder}/AMPsphere_GMSC_correspondence.tsv.gz')

    print('Analyzing families')
    spheres(analysis_folder)

    print('Generate fasta files')
    fasta_file(analysis_folder)
    
