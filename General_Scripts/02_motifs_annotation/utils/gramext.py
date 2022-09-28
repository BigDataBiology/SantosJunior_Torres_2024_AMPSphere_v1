import pandas as pd
import matplotlib.pyplot as plt
from itertools import product
from Bio import SeqIO


def convseq(seq):
    toA = ['A', 'D', 'E', 'K', 'R', 'N', 'T', 'S', 'Q']
    toC = ['C', 'F', 'L', 'I', 'V', 'M', 'W', 'Y', 'W', 'H']
    toT = ['P']
    for aa in toA: seq = seq.replace(aa, 'A')
    for aa in toC: seq = seq.replace(aa, 'C')
    for aa in toT: seq = seq.replace(aa, 'T')
    return seq


def gram_extract(seq):
    ngrams = [''.join(x) for x in product(['A','C','T','G'], repeat=4)]
    return [ seq.count(n) for n in ngrams ]


def bin_extract(seq):
    out = []
    ngrams = [''.join(x) for x in product(['A','C','T','G'], repeat=4)]
    for n in ngrams:
        if n in seq: out.append(1)
        else: out.append(0)
    return out
   

dfa, dfb = [], []
for record in SeqIO.parse('bacillus.fa', 'fasta'):
     record.seq = convseq(record.seq)
     sdf = gram_extract(record.seq)
     bdf = bin_extract(record.seq)
     sdf.insert(0, record.id)
     bdf.insert(0, record.id)
     dfa.append(sdf)
     dfb.append(bdf)

dfa = pd.DataFrame(dfa)
dfb = pd.DataFrame(dfb)

nms = [''.join(x) for x in product(['A','C','T','G'], repeat=4)]
nms.insert(0, 'seq')

dfa.columns = nms
dfb.columns = nms

dfa.to_csv('pseudo_code_bacillus.tsv.gz',
           sep='\t',
           header=True,
           index=None)

print(dfa)
print(dfb.drop('seq', axis=1).sum(axis=0).sort_values() / 70)

tfa, tfb = [], []
for record in SeqIO.parse('AMPSphere_v.2021-03.faa', 'fasta'):
     record.seq = convseq(record.seq)
     print(record.id)
     sdf = gram_extract(record.seq)
     bdf = bin_extract(record.seq)
     sdf.insert(0, record.id)
     bdf.insert(0, record.id)
     tfa.append(sdf)
     tfb.append(bdf)

tfa = pd.DataFrame(tfa)
tfb = pd.DataFrame(tfb)

tfa.columns = nms
tfb.columns = nms

tfa.to_csv('pseudo_code_AMPSphere.tsv.gz',
           sep='\t',
           header=True,
           index=None)

print(tfa)
print(tfb.drop('seq', axis=1).sum(axis=0).sort_values() / 863498)

i = dfb.drop('seq', axis=1).sum(axis=0).sort_values() / len(dfb)
ii = tfb.drop('seq', axis=1).sum(axis=0).sort_values() / len(tfb)

i.name, ii.name = 'bacillus', 'ampsphere'

df = pd. concat([i, ii], axis=1)
df.to_csv('trigram.tsv', sep='\t', header=True, index=True)

df.plot.bar()
plt.savefig('trigrams.png', dpi=300)

