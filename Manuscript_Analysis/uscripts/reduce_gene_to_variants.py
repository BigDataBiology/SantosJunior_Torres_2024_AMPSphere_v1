import lzma
import json
import pandas as pd

from Bio import SeqIO


def revseq(seq):
    stops = {'TTA', 'CTA', 'TCA'}
    if seq[0:3] in stops:
        return seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    else:
        return seq

 
hseqs = []
for r in SeqIO.parse(lzma.open('data/AMPSphere_v.2022-03.fna.xz',
                               'rt'),
                     'fasta'):
    access, amp = r.description.split(' | ')
    hseqs.append((access, amp, str(r.seq)))


hseqs = pd.DataFrame(hseqs, columns=['header', 'amp', 'seq'])
hseqs['seq'] = hseqs.seq.map(lambda x: revseq(x))
hseqs = hseqs.groupby(['seq', 'amp'])['header'].apply(lambda x: list(x))


hseqs = hseqs.reset_index()
hseqs = hseqs.sort_values(by=['amp', 'seq'])
hseqs = hseqs.reset_index(drop=True)
hseqs['amp_header'] = [x+'__'+str(j) for j, x in enumerate(hseqs.amp)]


redfasta = open('reduced_ampsphere_v2022-03.fna', 'wt')
for _, s, _, _, h in hseqs.itertuples(): redfasta.write(f'>{h}\n{s}\n')


rt = hseqs[['amp_header', 'header']]
final_dict = dict()
for _, h, gen in rt.itertuples():
    for i in gen: final_dict[i] = h


with open('ampsphere_gene_variant.json', 'w') as fp:
    json.dump(final_dict, fp)
