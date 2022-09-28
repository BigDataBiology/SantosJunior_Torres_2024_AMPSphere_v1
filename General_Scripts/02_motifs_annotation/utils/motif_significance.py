import re
import numpy as np
from Bio import SeqIO


def load_motif():
    import pandas as pd
    data = pd.read_table('data/db_motif.tsv',
                         sep='\t',
                         header='infer')
    dbmotif = data[['motif', 'name']]
    dbmotif = dbmotif.groupby('name')
    dbmotif = dbmotif.agg(lambda x: list(x))
    return dbmotif.to_dict()['motif']    
    
    
def anno_motif(sequence, dbmotif):
    motif_match = []
    for name, mlist in dbmotif.items():
        for motif in mlist:
            if re.search(motif, sequence):
                m = re.findall(motif, sequence)
                for x in m:
                    yield (name, x)


def pipe_anno(infile, ofile):
    dbmotif = load_motif()
    with open(ofile, 'w') as fout:
        fout.write('id\tmotif\tseq\n')
        for record in SeqIO.parse(infile, 'fasta'):
            record.seq = str(record.seq)
            for x in list(anno_motif(record.seq, dbmotif)):
                fout.write(f'{record.id}\t{x[0]}\t{x[1]}\n')
            print(f'{record.id}')
            
