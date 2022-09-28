def lfam_rnacode_output():
    import os
    import pandas as pd
    from glob import glob
    from collections import Counter
    from utils.genes_to_clusters import load_spheres

    address = 'analysis/large_clusters/rnacode/*.rnacode'

    df = pd.DataFrame()
    for i in glob(address):
        if os.stat(i).st_size != 0:
            f = i.split('/')[-1]
            f = f.split('_')[0:2]
            f = '_'.join(f)
            n = pd.read_table(i, header=None)
            n['family'] = f
            df = df.append(n)
    
    df.columns = ['HSS#', 'sense', 'frame',
                  'length', 'start_aa', 'end_aa', 
                  'rep_gene', 'start', 'end', 
                  'score', 'p', 'cluster']

    df = df[(df.sense == '+') & (df.frame == 1)]
    df.to_csv('analysis/large_clusters/rnacode_out.merged.tsv',
              sep='\t',
              header=True,
              index=None)

    selected_fams = Counter(df['cluster']).items()
    selected_fams = [k for k,v in selected_fams if v >= 2]

    fam = load_spheres()
    df = fam[fam.family.isin(selected_fams)]
    df = df.reset_index(drop=True)
    df['RNAcode_passed'] = 'True'
        
    flist = set()
    for i in glob(address):
        f = i.split('/')[-1].split('_')[0:2]
        f = '_'.join(f)
        flist.add(f)

    selected_fams = set(selected_fams)
    failed = flist - selected_fams
    failed = fam[fam.family.isin(failed)]
    failed = failed.reset_index(drop=True)
    failed['RNAcode_passed'] = 'False'

    df = df.append(failed)
    df = df.sort_values('amp')
    df.rename({'amp': 'AMP'}, axis=1, inplace=True)
    df.drop('family', axis=1, inplace=True)
    
    clean_names = df.AMP.tolist()
    prev_res = pd.read_table('analysis/RNAcode_amps_result.tsv')
    prev_res = prev_res[~prev_res.AMP.isin(clean_names)]
    prev_res = prev_res.append(df)
    prev_res = prev_res.sort_values('AMP')
    prev_res = prev_res.drop_duplicates()
    prev_res = prev_res.replace(['False', 'True', 'not_tested'],
                                ['Failed', 'Passed', 'Not tested'])

    prev_res.to_csv('RNAcode_out_wlfam.tsv', 
                    sep='\t', header=True,
                    index=None)
                    

