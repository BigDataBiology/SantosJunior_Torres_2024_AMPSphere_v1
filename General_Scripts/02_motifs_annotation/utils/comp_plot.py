def decomposition(seq, seqid):
    aas = 'ACDEFGHIKLMNPQRSTVYW'
    n = [seqid]
    for i in aas: n.append(seq.count(i))
    return n


def compplot():
    import gzip
    import pandas as pd
    import seaborn as sns
    from itertools import chain
    from collections import Counter
    from matplotlib import pyplot as plt
    from Bio import SeqIO
    from .environment_classification import environments_list

    envdict = dict()
    for ev, vv in zip(['non_mammalian', 'mammalian_gut', 'mammalian_other', 
                       'built_env.', 'freshwater', 'extreme',
                       'wastewater', 'soil', 'plant associated',
                       'marine', 'other'], environments_list):
        for i in vv: envdict[i] = ev
    

    genes = pd.read_table('data/gmsc_amp_genes_envohr_source.tsv.gz')
    genes = genes[genes.is_metagenomic == True]
    genes = genes[['amp', 'general_envo_name']]
    genes = genes.drop_duplicates()

    seqs = gzip.open('data/AMPSphere_v.2022-03.faa.gz',
                     'rt',
                     encoding='utf-8')

    ampsphere = []
    for record in SeqIO.parse(seqs, 'fasta'):
        ampsphere.append(decomposition(record.seq,
                                       record.id))

    columns = list('ACDEFGHIKLMNPQRSTVYW')
    columns.insert(0, 'amp')
    ampsphere = pd.DataFrame(ampsphere, columns=columns)
    
    genes = genes.merge(on='amp', right=ampsphere)
    genes = genes.fillna(0)
    genes.general_envo_name = [envdict[x] for x in genes.general_envo_name]
    genes.drop('amp', axis=1, inplace=True)

    df = genes.groupby('general_envo_name').agg('count')
    df = df['A']
    df = df.rename({'A': 'seqs'}, axis=1)
    
    genes = genes.groupby('general_envo_name').agg('sum')
    genes = genes.agg(lambda x: x/df)
    
    sns.clustermap(genes, cmap='YlOrBr', z_score=0)
    plt.savefig('analysis/amino_acids_composition.svg')
    plt.close()
    
    positive = genes[['R', 'H', 'K']].sum(axis=1)
    negative = genes[['D', 'E']].sum(axis=1)
    polar_uncharged = genes[['S', 'T', 'N', 'Q']].sum(axis=1)
    special_cases = genes[['G', 'P', 'C']].sum(axis=1)
    aromatic = genes[['F', 'W', 'Y']].sum(axis=1)
    hydrophobic = genes[['A', 'V', 'L', 'I', 'M']].sum(axis=1)
    grouped = pd.concat([positive,
                         negative,
                         polar_uncharged,
                         special_cases,
                         aromatic,
                         hydrophobic], axis=1)
    grouped.columns = ['positive',
                       'negative',
                       'polar_uncharged',
                       'special_cases',
                       'aromatic',
                       'hydrophobic']
    sns.clustermap(grouped, cmap='YlOrBr', z_score=0)
    plt.savefig('analysis/amino_acids_groupedI_composition.svg')
    plt.close()
    
    genes.to_csv('analysis/amino_acids_distribution_by_seq_by_env.tsv', sep='\t', header=True, index=None)
    grouped.to_csv('analysis/groupped_amino_acids_by_seq_by_env.tsv', sep='\t', header=True, index=None)
    


