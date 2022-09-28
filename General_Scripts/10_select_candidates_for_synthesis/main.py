def init_df():
    import pandas as pd

    print('load data')
    print('\tload features')
    features = pd.read_table('data/ampsphere_v2022-03.features.tsv.gz')
    features = features.rename({'Access': 'amp'}, axis=1).drop('group', axis=1).set_index('amp')
    print('\tload quality')
    quality = pd.read_table('data/quality_assessment.tsv')
    quality = quality.rename({'AMP': 'amp'}, axis=1).set_index('amp')
    print('\tload sol./syn. rules')
    solsyn = pd.read_table('data/AMPSphere_v.2022-03.solsyn_rules.tsv.gz')
    solsyn = solsyn.rename({'accession': 'amp'}, axis=1).drop('fam', axis=1).set_index('amp')
    print('\tload families')
    families = pd.read_table('data/SPHERE_v.2022-03.levels_assessment.tsv.gz')
    families = families.rename({'AMP accession': 'amp', 'SPHERE_fam level III': 'family'}, axis=1)[['amp', 'family']].set_index('amp')
    print('\tload coprediction')
    copred = pd.read_table("data/ensemble_predictions.tsv")
    print('\tload gmsc_ampsphere data')
    data = pd.read_table("data/gmsc_amp_genes_envohr_source.tsv.gz")

    print('organize data')
    print('\tselect data info')
    copies = data.groupby('amp').agg('size')
    d2 = data[['amp', 'source']].drop_duplicates().fillna('*').groupby('amp')['source'].apply(lambda x: ', '.join(set(x)))
    envo = data[['amp', 'general_envo_name']].fillna('isolated genome').groupby('amp')['general_envo_name'].apply(lambda x: ', '.join(set(x)))
    samp = data[['amp', 'sample']].drop_duplicates().groupby('amp')['sample'].apply(lambda x: ', '.join(set(x)))
    newdf = pd.concat([copies, d2, envo, samp], axis=1).rename({0: 'copies'}, axis=1)

    print('\tselect features info')
    newdf = pd.concat([newdf, features], axis=1)
    newdf = newdf.drop(['SA.Group1.residue0', 'SA.Group2.residue0',
                    'SA.Group3.residue0', 'HB.Group1.residue0',
                    'HB.Group2.residue0', 'HB.Group3.residue0'], axis=1)
    newdf = newdf.drop(['tinyAA', 'smallAA', 'aliphaticAA',
                    'aromaticAA', 'nonpolarAA', 'polarAA',
                    'chargedAA', 'basicAA', 'acidicAA'], axis=1)

    print('\tselect quality and families info')
    newdf = pd.concat([newdf, quality], axis=1)
    newdf = pd.concat([newdf, families], axis=1)

    print('\tselect fix info')
    newdf['length_residues'] = newdf.sequence.apply(lambda x: len(x))
    newdf = newdf.rename({'Coordinates': 'terminal placement'}, axis=1)

    newdf = newdf[['family', 'copies', 'source',
                   'general_envo_name', 'sample', 'sequence',
                   'length_residues', 'charge', 'pI',
                   'aindex', 'instaindex',
                   'boman', 'hydrophobicity',
                   'hmoment', 'Antifam',
                   'RNAcode', 'metaproteomes',
                   'metatranscriptomes', 'terminal placement']]

    print('\tadd final data')
    newdf = pd.concat([newdf, solsyn], axis=1)
    newdf = newdf.reset_index().rename({'general_envo_name': 'environment'}, axis=1)
    return newdf, copred

    
def output_df(newdf, copred):
    print('output data')
    print('\tsaving sequences')
    newdf[['amp', 'family',
           'length_residues',
           'sequence']].to_csv('analysis/AMPSphere_sequences.tsv',
                               sep='\t',
                               header=True,
                               index=None)
    print('\tsaving taxonomy')
    newdf[['amp', 'copies', 'source',
           'environment', 'sample']].to_csv('analysis/AMPSphere_taxonomy_environment.tsv',
                                            sep='\t',
                                            header=True,
                                            index=None)
    print('\tsaving features')
    newdf[['amp', 'charge', 'pI',
           'aindex', 'instaindex',
           'boman', 'hydrophobicity',
           'hmoment']].round(3).to_csv('analysis/AMPSphere_features.tsv',
                                       sep='\t',
                                       header=True,
                                       index=None)
    print('\tsaving qc')
    newdf[['amp', 'Antifam', 'RNAcode',
           'metaproteomes', 'metatranscriptomes',
           'terminal placement']].to_csv('analysis/AMPSphere_quality_control.tsv',
                                         sep='\t',
                                         header=True,
                                         index=None)
    print('\tsaving sunsol')
    newdf[['amp', 'sol_rule_1', 'sol_rule_2',
           'sol_rule_3', 'sol_rule_4', 'sol_rule_5',
           'sol_rule_6', 'syn_rule_1', 'syn_rule_2',
           'syn_rule_3']].to_csv('analysis/AMPSphere_synthesis_solubility.tsv',
                                 sep='\t',
                                 header=True,
                                 index=None)
    print('\tsave copred')
    copred.to_csv('analysis/AMPSphere_selected_candidates_coprediction.tsv',
                  sep='\t',
                  header=True,
                  index=None)


def co_loading():
    import pandas as pd
    from itertools import chain
    # loading homologs
    smprot = pd.read_table('data/result_SmProt.m8',
                           sep='\t',
                           header=None)
    starpep = pd.read_table('data/result_starPepDB.m8',
                            sep='\t',
                            header=None)
    dramp = pd.read_table('data/result_dramp.m8',
                          sep='\t',
                          header=None)
    homologs = pd.concat([dramp, starpep, smprot])
    homologs = homologs[homologs[2] <= 1e-5]
    homologs.columns = ['query', 'target', 'evalue',
                        'gapopen', 'pident', 'nident',
                        'qstart', 'qend', 'qlen',
                        'tstart', 'tend', 'tlen',
                        'alnlen', 'raw', 'bits',
                        'cigar', 'qseq', 'tseq',
                        'qheader', 'theader', 'qaln',
                        'taln', 'qframe', 'tframe',
                        'mismatch', 'qcov', 'tcov']
    homologs.to_csv('analysis/AMPSphere_homologs.tsv',
                    sep='\t',
                    header=True,
                    index=None)
    # getting homologs to all tested dbs
    homologs = set(homologs['query'])
    # load fmt
    fmt = pd.read_table("data/AMPs_per_genomes.tsv")
    fmt = fmt[~fmt.sequences.isna()]
    fmt.sequences = fmt.sequences.apply(lambda x: x.split(","))
    fmt_amps = set(chain.from_iterable(fmt.sequences))
    return homologs, fmt_amps


def gen_filter(newdf, copred):
    from scipy.stats.mstats import gmean
    # selection of first batch, 50 peps
    # 1. max len -- 40 res.
    l = (newdf.length_residues <= 40)
    # 2. being approved in at least 60% of sol/syn tests
    newdf['sol']= newdf[['sol_rule_1',
                         'sol_rule_2',
                         'sol_rule_3',
                         'sol_rule_4',
                         'sol_rule_5',
                         'sol_rule_6']].sum(axis=1) / 6
    newdf['syn'] = newdf[['syn_rule_1',
                          'syn_rule_2',
                          'syn_rule_3']].sum(axis=1) / 3
    newdf['fscore'] = gmean(newdf[['sol', 'syn']], axis=1)
    sc = (newdf.fscore >= 0.6)
    # pre-selecting
    df = newdf[l & sc]
    # 3. being co-predicted by all methods
    df = df[df.amp.isin(copred[copred.pcts == 100]['amp'])]
    red_df = df.sort_values(by=['family',
                                'length_residues',
                                'copies'],
                            ascending=[True,
                                       False,
                                       False])
    red_df = red_df.groupby('family').apply(lambda x: x.head(1))
    return red_df

def main():
    newdf, copred = init_df()
    out_df(newdf, copred)
    red_df = gen_filter(newdf, copred)
    df = red_df.sort_values(by=['copies',
                                'fscore',
                                'length_residues'],
                            ascending=[False,
                                       False,
                                       False])
    df = df[['amp', 'family', 'copies', 'fscore', 'sequence']]
    df.reset_index(drop=True, inplace=True)
    homologs, fmt_amps = co_loading()
    df['feature'] = ''
    df.loc[df.copies >= 50, 'feature'] = 'high-copy'
    df.loc[df.copies < 50, 'feature'] = 'low-copy'
    df.loc[df.sequence.isin(homologs), 'feature'] = 'homolog'
    df.loc[df.amp.isin(fmt_amps), 'feature'] = 'fmt-related'
    df = df.merge(on='amp', right=newdf[['amp', 'source', 'environment']])
    df.to_csv('analysis/selected_candidates.tsv',
              sep='\t',
              header=True,
              index=None)

if __name__ == '__main__':
    main()
    
