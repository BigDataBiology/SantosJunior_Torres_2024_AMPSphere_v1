import pandas as pd
import numpy as np
from jug import TaskGenerator


@TaskGenerator
def recalculated_evalue(level, gap_open, gap_extend):
    import pandas as pd
    from utils import clustering_analysis
    from Bio.pairwise2 import align
    import Bio.Align.substitution_matrices as mt
    import numpy as np
    mapped = pd.read_table(f'../../data_folder/output_clustering_significance_level{level}.tsv.gz')
    sequences = clustering_analysis.seqload()

    re_calc_evalue = []
    for ix in range(len(mapped)):
        seq1 = sequences[mapped.iloc[ix]['query']]
        seq2 = sequences[mapped.iloc[ix]['target']]

        if gap_open < 0 or gap_extend < 0:
            alignments = align.globalds(seq1,
                                        seq2,
                                        mt.load('BLOSUM62'),
                                        gap_open,
                                        gap_extend)
        else:
            alignments = align.globaldx(seq1,
                                        seq2,
                                        mt.load('BLOSUM62'))
        aligned_A, aligned_B, score, begin, end = alignments[0]
        re_calc_evalue.append(clustering_analysis.f_evalue(seq1, score, (gap_open != 0 or gap_extend != 0)))

    return np.array(re_calc_evalue)

@TaskGenerator
def summarize(r):
    return {k: np.mean(v < 1e-5) for k,v in r.items()}

@TaskGenerator
def write_out(final):
    from os import makedirs
    makedirs('outputs', exist_ok=True)
    pd.Series(final).reset_index().to_csv('outputs/gap-penalty-comparison.tsv', sep='\t', index=False, header=False)

results = {}
for level in ['I', 'II', 'III']:
    for gap_open, gap_extend in [
        (-10, -.5),
        ( -9,  -1),
        (-11,  -1),
        ( -5, -.5),
        ( -5,  -1),
        ( -1,  -1),
        (  0,   0),
        ]:
        results[level, gap_open, gap_extend] = recalculated_evalue(level, gap_open, gap_extend)

final = summarize(results)
write_out(final)
