def main():
    import os
    import gzip
    import pandas as pd
    
    from Bio import SeqIO
    from scipy.stats.mstats import gmean
    from utils.download_files import inputsgen
    from utils.sol_synth import solubility, synthesis

    print('# generate folders')
    os.makedirs('data/', exist_ok=True)
    os.makedirs('analysis/', exist_ok=True)

    print('# download files')
    inputsgen()
    
    print('# stating files')
    fin = gzip.open('data/AMPSphere_v.2022-03.faa.gz',
                    'rt', encoding='utf-8')
                    
    fam = 'data/SPHERE_v.2022-03.levels_assessment.tsv.gz'
    fam = pd.read_table(fam, sep='\t', header='infer')
    fam.rename({'AMP accession': 'accession',
                'SPHERE_fam level III': 'fam'},
               axis=1, inplace=True)
    fam = fam[['accession', 'fam']] 
              
    hq = 'data/high_quality_candidates.txt'
    hq = pd.read_table(hq, header=None)
    
    q = 'data/quality_candidates.txt'
    q = pd.read_table(q, header=None)
    
    quality_candidates = pd.concat([hq, q])[0]
    
    print('# checking solubility and synthesis rules')
    res = []
    for record in SeqIO.parse(fin,
                              'fasta'):
        sol_val = solubility(str(record.seq))
        sol_val = sol_val.rule_judgement()
        syn_val = synthesis(str(record.seq))
        syn_val = syn_val.rule_judgement()
        res.append([record.id] + sol_val + syn_val)

    print('# converting in DF')        
    evaluation = pd.DataFrame(res, columns=['accession',
                                            'sol_rule_1',
                                            'sol_rule_2',
                                            'sol_rule_3',
                                            'sol_rule_4',
                                            'sol_rule_5',
                                            'sol_rule_6',
                                            'syn_rule_1',
                                            'syn_rule_2',
                                            'syn_rule_3'])
    
    print('# adding up families')
    evaluation = evaluation.merge(on='accession',
                                  right=fam)
    
    evaluation.to_csv('analysis/AMPSphere_v.2022-03.solsyn_rules.tsv.gz',
                      sep='\t', header=True, index=None)
    
    print('# adding up quality')
    evaluation = evaluation[evaluation.accession.isin(quality_candidates)]

    evaluation.to_csv('analysis/AMPSphere_v.2022-03.selected_candidates.tsv.gz',
                      sep='\t', header=True, index=None)                      
    
    print('# calculating scores for the best candidates')
    evaluation['sol_score'] =   evaluation[['sol_rule_1',
                                            'sol_rule_2',
                                            'sol_rule_3',
                                            'sol_rule_4',
                                            'sol_rule_5',
                                            'sol_rule_6']].T.sum() * 100 / 6

    evaluation['syn_score'] =   evaluation[['syn_rule_1',
                                            'syn_rule_2',
                                            'syn_rule_3']].T.sum() * 100 / 3

    evaluation['FScore'] = gmean(evaluation[['sol_score',
                                             'syn_score']],
                                 axis=1)
                                 
    evaluation = evaluation[evaluation.FScore >= 60]
    evaluation = evaluation.drop(['sol_score',
                                  'syn_score'],
                                 axis=1)
    
    evaluation = evaluation.sort_values(by=['FScore',
                                            'accession',
                                            'fam'],
                                        ascending=[False,
                                                   True,
                                                   True])

    evaluation.to_csv('analysis/selected_candidates.tsv',
                      sep='\t', header=True, index=None)
                      
                      
if __name__ == '__main__':
    main()

