def load_aln():
    import pandas as pd
    df = pd.read_table('analysis/file_order.tsv', 
                       sep='\t', header='infer')
    return df
    
    
def merge_rnacode():
    import os
    import pandas as pd
    from glob import glob
    df = pd.DataFrame()
    for i in glob('analysis/rnacode_out/*rnacode'):
        if os.stat(i).st_size != 0:
            f = i.split('/')[-1].replace('.rnacode', '')
            n = pd.read_table(i, header=None)
            n['cluster'] = f
            df = pd.concat([df, n])
    
    df.columns = ['HSS#', 'sense', 'frame',
                  'length', 'start_aa', 'end_aa', 
                  'rep_gene', 'start', 'end', 
                  'score', 'p', 'cluster']
    
    df = df[(df.sense == '+') & (df.frame == 1)]
    aln = load_aln()
    df = df.merge(on='cluster', right=aln)
    df.to_csv('analysis/rnacode_output.tsv',
              sep='\t',
              header=True,
              index=None)
    

def sumlist(subset):
     n = []
     for i in subset:
         i = i.split(',')
         n += i
     return n
     

def list_results():
    import pandas as pd
    from collections import Counter
    from .genes_to_clusters import load_spheres
    true_list = pd.read_table('analysis/rnacode_output.tsv',
                              sep='\t', header='infer')
    spotted_fams = true_list.cluster.tolist()
    true_list = sumlist(true_list.amps)
    aln = load_aln()
    aln = aln[~aln.cluster.isin(spotted_fams)]
    false_list = sumlist(aln.amps)
    clusters = load_spheres()
    f1 = (~clusters.amp.isin(true_list))
    f2 = (~clusters.amp.isin(false_list))
    not_tested = clusters[f1 & f2]['amp']
    out = []
    for i in true_list: out.append([i, 'Passed'])
    for i in false_list: out.append([i, 'Failed'])
    for i in not_tested: out.append([i, 'Not tested'])
    out = pd.DataFrame(out, columns=['AMP', 'RNAcode_passed'])
    out = out.sort_values(by='AMP')
    print(Counter(out.RNAcode_passed))
    out.to_csv('analysis/RNAcode_amps_result.tsv',
               sep='\t', header=True, index=None)
               
