def get_orfs(folder):
    '''
    Get the total number of ORFs and smORFs
    Retrieve the info about 
    '''
    import pandas as pd
    from Bio import SeqIO
    from os.path import exists        
    # getting amps
    f = f'{folder}/macrel.out.prediction.gz'
    if exists(f):
       data = pd.read_table(f,
                            sep='\t',
                            header='infer',
                            comment='#')
       amps = set(data['Access'])
       seqs = ','.join(set(data['Sequence']))
    else: amps, seqs = {}, ''
    n_amps = len(amps)
    # getting orfs                        
    smorfs, orfs = [0, 0]
    f = f'{folder}/macrel.out.all_orfs.faa'
    for record in SeqIO.parse(f, 'fasta'):
        orfs += 1
        if len(record.seq) <= 100:
            smorfs += 1
    return orfs, smorfs, n_amps, seqs
            
            
def eval_results():
    import glob
    import pandas as pd
    metadata = pd.read_table('metadata/selected_genomes.tsv')
    df = []
    for idx, folder in enumerate(glob.glob('macrel_results/*_macrel')):
        genome = folder.split('/')[-1].replace('_macrel', '')
        orfs, smorfs, n_amps, seqs = get_orfs(folder)
        df.append([genome, orfs, smorfs, n_amps, seqs])
        print(idx)
    df = pd.DataFrame(df,
                      columns=['genome',
                               'ORFs',
                               'smORFs',
                               '#_amps', 
                               'sequences'])
    df = df.merge(on='genome',
                  right=metadata,
                  how='outer')
    df.to_csv('AMPs_per_genomes.tsv',
              sep='\t',
              header=True,
              index=None)
    df = df[['ORFs', 'smORFs',
             '#_amps', 'cluster',
             'subject_type',
             'timepoint.fmt',
             'clinical_response']]
    df = df.fillna('n.a.')
    df = df.pivot_table(index=['cluster',
                               'subject_type',
                               'timepoint.fmt',
                               'clinical_response'],
                        aggfunc=(lambda x: list(x)))
    df = df.reset_index()
    df = df.to_csv('amps_cluster_out.tsv', sep='\t', header=True, index=None)

