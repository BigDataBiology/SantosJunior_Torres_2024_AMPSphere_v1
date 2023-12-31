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
        seqdict[record.id] = str(record.seq)
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
        for _, record in d.groupby(f'SPHERE_fam level {l}'):
            record = record[record.L == record.L.max()]
            if len(record) > 1:
                record = record.sort_values(by='sequence')
            o.append(record['AMP accession'].iloc[0])
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
    return tobesampled.sample(n)['AMP accession'].tolist()


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


def aln(seq1,
        seq2,
        gap_open=-10,
        gap_extend=-0.5):
    '''
    Align two sequences using BLOSUM62
    Returns score and alignment length
    '''
    from Bio.pairwise2 import align
    import Bio.Align.substitution_matrices as mt
    alignments = align.globalds(seq1,
                                seq2,
                                mt.load('BLOSUM62'),
                                gap_open,
                                gap_extend)
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
    K, l = (0.041, 0.267)
    # Number of residues in the AMPSphere dataset
    m = 32509722
    N = m*len(seq1)
    Sbit = l*score
    Sbit -= np.log(K)
    Sbit /= np.log(2)
    # The score is an HSP, then the evalue formula is: Kmne^(-l*S)
    # The bit score (Sbit) use is more accurate because it
    # was optimized. The e-value is calculated as: N/(e^Sbit)
    return N * np.exp(-Sbit)


def align_to_representative(sample_df):
    '''
    Process the samples by replicate returning
    a dataframe of the alignment of each sampled sequence
    against their cluster representative.
    '''
    import pandas as pd
    out = []
    for record in sample_df.itertuples():
        sid, gid, cov, score, lent = aln(record.sequence_x,
                                         record.sequence_y, -10, -0.5)
        evalue = f_evalue(record.sequence_x, score)
        out.append([record._1,
                    record._4,
                    sid,
                    gid,
                    cov,
                    score,
                    lent,
                    evalue,
                    record.family])
    out = pd.DataFrame(out, columns=['query',
                                     'target',
                                     'identity',
                                     'gap_identity',
                                     'coverage',
                                     'score',
                                     'aln_len',
                                     'evalue',
                                     'family'])
    out['sig.'] = out['evalue'].apply(lambda x: '*' if x < 1e-5 else 'n.s.')
    return out


def cluster_analysis():
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

        sampled = sample_seqs(clusters,
                                 rep,
                                 level,
                                 3000)
        sampled = clusters[clusters['AMP accession'].isin(sampled)]
        cols = ['AMP accession',
                  f'SPHERE_fam level {level}',
                  'sequence']
        sampled = sampled[cols]
        repdata = clusters[clusters['AMP accession'].isin(rep)]
        sampled = sampled.merge(on=f'SPHERE_fam level {level}',
                       right=repdata[cols])
        sampled.rename(columns={f'SPHERE_fam level {level}': 'family'},
                        inplace=True)
        sampled = align_to_representative(sampled)
        sampled.to_csv(f'output_clustering_significance_level{level}.tsv',
                  sep='\t',
                  header=True,
                  index=None)

