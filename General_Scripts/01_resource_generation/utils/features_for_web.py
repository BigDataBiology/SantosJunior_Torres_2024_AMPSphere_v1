def protein_check(seq):
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    x = ProteinAnalysis(seq)
    mw = x.molecular_weight()
    ar = x.aromaticity()
    ii = x.instability_index()
    pi = x.isoelectric_point()
    charge = x.charge_at_pH(7.0)
    gravy = x.gravy()
    return [seq, mw, ar, ii, pi, charge, gravy]


def features_for_web():
    import gzip
    import pandas as pd
    from Bio import SeqIO

    fams = pd.read_table('analysis/SPHERE_v.2022-03.levels_assessment.tsv.gz')
    fams = fams[['AMP accession', 'SPHERE_fam level III']]
    fams.columns = ['id', 'family']

    seqs = []
    for record in SeqIO.parse(gzip.open('analysis/AMPSphere_v.2022-03.faa.gz',
                                        'rt'),
                              'fasta'):
        seqs.append((record.id, record.seq))

    seqs = pd.DataFrame(seqs, columns=['id', 'sequence'])
    seqs.sequence = [''.join(x) for x in seqs.sequence]
    seqs['length'] = seqs.sequence.apply(lambda x: len(x))

    prot_feat = [protein_check(x) for x in seqs.sequence]

    prot_feat = pd.DataFrame(prot_feat, columns=['sequence',
                                                 'molecular_weight',
                                                 'aromaticity',
                                                 'instability_index',    
                                                 'isoelectric_point',
                                                 'charge',
                                                 'gravy'])

    seqs = seqs.merge(on='sequence', right=prot_feat)
    seqs = seqs.drop('sequence', axis=1)
    seqs = seqs.merge(on='id', right=fams)
    seqs = seqs[['id',
                 'family',
                 'length',
                 'molecular_weight',    
                 'aromaticity',
                 'instability_index',
                 'isoelectric_point',
                 'charge',    
                 'gravy']]

    seqs.to_csv('analysis/AMPSphere_v.2022-03.features_for_web.tsv.gz',
                sep='\t',
                header=True,
                index=None)    


def main():
    features_for_web()


if __name__=='__main__':
    main()
    
