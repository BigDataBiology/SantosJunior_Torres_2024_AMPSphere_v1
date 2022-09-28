import sys
import numpy as np
import pandas as pd

from Bio import SeqIO
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint as IP


def serum_half_time(sequence):
    '''
    ## This method of regression
    ## adds up to the half life estimation
    ## of peptides in the serum
    DOI: 10.1111/cts.12985
    '''
    sequence = sequence.upper()

    NP = 0
    for aa in ['A', 'L', 'I', 'P', 'W', 'V', 'F', 'M']:
        NP += sequence.count(aa)
    NP = NP * 100/ len(sequence)

    if sequence.count('W') > 0: w = 1
    else: w = 0

    if sequence.count('Y') < 2: y = 0
    else: y = 1

    if IP(sequence).pi() < 10: pI = 0
    else: pI = 1

    ln_hlf_t = 2.226 + 0.053*NP - 1.515*w + 1.290*y - 1.052*pI

    return round(np.exp(ln_hlf_t), 3)


def main(infile):
    headers, hl = [], []
    
    for record in SeqIO.parse(infile, 'fasta'):
        headers.append(record.id)
        hl.append(serum_half_time(record.seq))
        print(record.id)
        
    df = pd.DataFrame(np.array([headers, hl]).T)
    df.columns = ['sequence', 'half_life']
    df.to_csv('serum_half_lifes.tsv', sep='\t', header=True, index=None)


if __name__ == '__main__':
    main(sys.argv[1])
    
