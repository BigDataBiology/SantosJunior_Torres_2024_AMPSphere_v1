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


def seqload():
    '''
    Load AMPSphere peptide sequences
    as a dictionary
    '''
    import gzip
    from Bio import SeqIO
    f = gzip.open('data/AMPSphere_v.2022-03.faa.gz',
                  'rt', encoding='utf-8')
    seqdict = dict()
    for record in SeqIO.parse(f, 'fasta'):
        record.seq = reduceseq(str(record.seq))
        seqdict[record.id] = record.seq
    return seqdict


def clusters_load():
    '''
    Loads the cluster info as a dataframe
    adding up sequence and length to the dataset
    '''
    import pandas as pd
    # load files
    data = pd.read_table('data/SPHERE_v.2022-03.levels_assessment.tsv.gz')
    seqdict = seqload()
    # get sequence data
    data['sequence'] = [seqdict[x] for x in data['AMP accession']]
    data['L'] = [len(x) for x in data.sequence]
    return data


def rep_level(data):
    '''
    Select the representative sequences for each non-singleton
    cluster at the three different clustering levels
    '''
    # select representatives
    from collections import Counter
    rI, rII, rIII = [], [], []
    dl_list = [['I', rI],
               ['II', rII],
               ['III', rIII]]
    for l, o in dl_list:
        print(f'... level {l}')
        excl = data[f'SPHERE_fam level {l}']
        excl = Counter(excl).items()
        excl = [k for k, v in excl if v > 1]
        filt = ['AMP accession',
                f'SPHERE_fam level {l}',
                'sequence',
                'L']
        d = data[filt]
        d = d[d[f'SPHERE_fam level {l}'].isin(excl)]
        for record in d.groupby(f'SPHERE_fam level {l}'):
            record = record[1]
            record = record[record.L == record.L.max()]
            if len(record) == 1:
                o.append(record['AMP accession'].tolist()[0])
            else:
                record = record.sort_values(by='sequence')
                record = record.head(1)
                o.append(record['AMP accession'].tolist()[0])
    return (rI, rII, rIII)


def sample_seqs(data, representatives: list, level: str, n=1000):    
    '''
    Draw samples for the level of clustering studied
    Inputs: data - dataframe obtained from clusters_load()
            representatives - list of AMP accessions of representatives
            level - string representing the clustering level (I, II, III)
    '''
    from collections import Counter
    # filter representatives off
    kexcl = Counter(data[f'SPHERE_fam level {level}']).items()
    kexcl = [k for k, v in kexcl if v > 1]
    tobesampled = data[(~data['AMP accession'].isin(representatives))]
    tobesampled = tobesampled[tobesampled[f'SPHERE_fam level {level}'].isin(kexcl)]
    tobesampled = tobesampled[['AMP accession', f'SPHERE_fam level {level}']]
    # perform sampling
    sample_1 = tobesampled.sample(n)['AMP accession'].tolist()
    sample_2 = tobesampled[(~tobesampled['AMP accession'].isin(sample_1))]
    sample_2 = sample_2.sample(n)['AMP accession'].tolist()
    sample_3 = tobesampled[(~tobesampled['AMP accession'].isin(sample_1 + sample_2))]
    sample_3 = sample_3.sample(n)['AMP accession'].tolist()
    return [sample_1, sample_2, sample_3]


def fix_rnsamples(data, sample, representatives, level):
    '''
    Create a final file to input the alignment test
    from representatives and samples from the sample
    clustering level
    Inputs: data is a dataframe from the clusters_load function
            sample_x are lists of accessions sampled for the test
            representatives is accession list of representatives
            for that level
            level is the clustering level, a string (I, II, III)
    '''
    cols = ['AMP accession',
            f'SPHERE_fam level {level}',
            'sequence']
    representatives = data[data['AMP accession'].isin(representatives)]
    sample = data[data['AMP accession'].isin(sample)]
    representatives = representatives[cols]
    sample = sample[cols]
    sample = sample.merge(on=f'SPHERE_fam level {level}',
                          right=representatives)
    sample.rename({f'SPHERE_fam level {level}': 'family'},
                  axis=1,
                  inplace=True)
    return sample


def _calculate_identity(sequenceA, sequenceB):
    """
    Returns the percentage of identical characters between two sequences.
    Assumes the sequences are aligned.
    """
    sa, sb, sl = sequenceA, sequenceB, len(sequenceA)
    matches = [sa[i] == sb[i] for i in range(sl)]
    length = len([1 for i in range(sl) if sa[i] != '-'])
    seq_id = (100 * sum(matches)) / sl
    gapless_sl = sum([1 for i in range(sl) if (sa[i] != '-' and sb[i] != '-')])
    gap_id = (100 * sum(matches)) / gapless_sl
    cov = (100 * gapless_sl) / length
    return (seq_id, gap_id, cov)


def aln(seq1, seq2):
    '''
    Align two sequences using BLOSUM62
    Returns score and alignment length
    '''
    from Bio.pairwise2 import align
    import Bio.Align.substitution_matrices as mt
    alignments = align.globalds(seq1,
                                seq2,
                                mt.load('BLOSUM62'),
                                -10,
                                -0.5)
    aligned_A, aligned_B, score, begin, end = alignments[0]
    lent = end - begin  # alignment length
    seq_id, g_seq_id, cov = _calculate_identity(aligned_A, aligned_B)
    return (seq_id, g_seq_id, cov, score, lent)


def f_evalue(seq1, score):
    '''
    From a sequence length and a score of the 
    alignment, calculates the corresponding evalue
    for hit
    '''
    import math
    import numpy as np
    # From BLAST adapted to peptides search (Lambda)
    K, l = (0.132539, 0.313667)
    # Number of residues in the AMPSphere dataset
    m = 32509722
    # start calculus
    N = m*len(seq1)
    Sbit = l*score
    Sbit = Sbit - np.log(K)
    Sbit = Sbit / np.log(2)
    # The score is an HSP, then the evalue formula is: Kmne^(-l*S)
    # The bit score (Sbit) use is more accuracte because it
    # was optimized. The e-value is calculated as: N/(e^Sbit)
    evalue = -1 * Sbit
    evalue = (math.e) ** evalue
    evalue = N * evalue
    if evalue <= 1e-5: ev = '*'
    else: ev = 'n.s.'
    return (evalue, ev)


def process_aln(sample_df, replicate):
    '''
    Process the samples by replicate returning
    a dataframe of the alignment of each sampled sequence
    against their cluster representative.
    '''
    import pandas as pd
    out = []
    for record in sample_df.itertuples():
        sid, gid, cov, score, lent = aln(record.sequence_x,
                                         record.sequence_y)
        evalue, ev = f_evalue(record.sequence_x,
                              score)
        out.append([record._1, 
                    record._4, 
                    sid,
                    gid, 
                    cov,
                    score, 
                    lent, 
                    evalue,
                    ev,
                    record.family])
    out = pd.DataFrame(out, columns=['query',
                                     'target',
                                     'identity',
                                     'gap_identity',
                                     'coverage',
                                     'score',
                                     'aln_len',
                                     'evalue',
                                     'sig.',
                                     'family'])
    out['replicate'] = replicate
    return out
    

def cluster_redalph():
    import pandas as pd
    print('load info')
    clusters = clusters_load()
    print('select representatives')
    rI, rII, rIII = rep_level(clusters)
    print('process by level')
    for rep, level in [[rI, 'I'],
                       [rII, 'II'],
                       [rIII, 'III']]:
        print(f'- sampling at level {level}')
        s1, s2, s3 = sample_seqs(clusters,
                                 rep,
                                 level,
                                 1000)
        print(f'... aligning first replicate')
        s1 = fix_rnsamples(clusters,
                           s1,
                           rep,
                           level)
        s1 = process_aln(s1,
                         '1')    
        print(f'.... aligning second replicate')
        s2 = fix_rnsamples(clusters,
                           s2,
                           rep,
                           level)
        s2 = process_aln(s2,
                         '2')    
        print(f'... aligning third replicate')
        s3 = fix_rnsamples(clusters,
                           s3,
                           rep,
                           level)
        s3 = process_aln(s3,
                         '3')
        print('... merging results')
        df = pd.concat([s1, s2, s3])
        print('... exporting')
        df.to_csv(f'output_clustering_redalph_significance_level{level}.tsv',
                  sep='\t',
                  header=True,
                  index=None)
    
