from macrel.fasta import fasta_iter
from macrel.AMP_features import fasta_features

with open('outputs/NAMPs.fasta', 'wt') as ofile:
    for h,seq in fasta_iter('AMP.train.fasta'):
        name,group = h.split('|')
        if group == 'NAMP':
            ofile.write(f'>{name}\n{seq}\n')

feat = fasta_features('outputs/NAMPs.fasta')
feat = feat.reset_index().rename(columns={'index': 'Access'})
feat.to_csv('new_data/macrel_trainneg.features.tsv.gz', sep='\t', header=True, index=None)
