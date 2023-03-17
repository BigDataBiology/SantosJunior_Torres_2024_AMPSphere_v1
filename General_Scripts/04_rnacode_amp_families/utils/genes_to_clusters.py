def revcomp(seq):
    nts = {'A': 'T', 'T': 'A',
           'C': 'G', 'G': 'C'}
    # reverting
    seq = seq[::-1]
    # complementing
    seq = [nts[x] for x in seq]
    # uniting
    seq = ''.join(seq)
    return seq
    
    
def check_sense(seq):
    stop_codons = ['TAA', 'TGA', 'TAG']
    terminal = seq[-3:]    
    if terminal in stop_codons: return False
    else: return True
    
          
def load_nts():
    import lzma
    import pandas as pd
    from Bio import SeqIO

    infile = lzma.open('data/AMPSphere_v.2022-03.fna.xz',
                       'rt',
                       encoding='utf-8')
    data = []
    for record in SeqIO.parse(infile, 'fasta'):
        print(record.description)
        record.seq = str(record.seq)
        n = record.description.split(' | ')
        if check_sense(record.seq):
            record.seq = revcomp(record.seq)
        data.append([n[0],
                     n[1],
                     str(record.seq)])

    data = pd.DataFrame(data,
                        columns=['gmsc',
                                 'amp',
                                 'sequence']) 
    return data
    

def load_spheres():
    import pandas as pd
    df = pd.read_table('data/SPHERE_v.2022-03.levels_assessment.tsv.gz')
    df = df[['AMP accession', 'SPHERE_fam level III']]
    df.rename({'AMP accession': 'amp', 'SPHERE_fam level III': 'family'},
              axis=1, inplace=True)
    return df
    
    
def merge_clusters():
    from collections import Counter
    genes = load_nts()
    clusters = load_spheres()
    genes = genes.merge(on='amp',
                        right=clusters)
    genes = genes[['amp',
                   'family',
                   'sequence']]
    genes = genes.drop_duplicates()
    genes['names'] = [f'{x}_{str(j)}' for j, x in enumerate(genes['amp'])]
    selected_fams = Counter(genes.family).items()
    selected_fams = [k for k,v in selected_fams if 500 >= v >= 3 ]
    genes = genes[genes.family.isin(selected_fams)][['family', 'names', 'sequence']]
    return genes


def align(in_file):
    import os
    from Bio import AlignIO
    from Bio.Align.Applications import MafftCommandline
    mafft_cline = MafftCommandline(input=in_file)
    stdout, stderr = mafft_cline()
    o_file = in_file.replace('fasta', 'aln')
    with open('o.aln', 'w') as odb:
        odb.write(stdout)
    AlignIO.convert('o.aln', 'fasta',
                    o_file, 'clustal')
    for f in ['o.aln', in_file]: os.remove(f)
        

def create_fasta(genes):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO
    for record in genes.groupby('family'):
        print(record[0])
        seqs = []
        for (idx, fam, name, seq) in record[1].itertuples():
            seqs.append(SeqRecord(Seq(seq),
                                  id=name,
                                  annotations={'molecule_type': 'DNA',
                                               'family': record[0]}))
        SeqIO.write(seqs, f'{record[0]}.fasta', 'fasta')
        align(f'{record[0]}.fasta')

