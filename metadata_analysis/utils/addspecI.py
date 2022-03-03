def addspecI():
    '''
    Add specI cluster info in the complete table of AMP, genes, and taxonomy
    '''
    import lzma, gzip
    import pandas as pd
    import numpy as np
    from tqdm import tqdm

    # define taxonomy levels
    tab_levels = {'genus': 'G',
                  'order': 'O',
                  'species': 'S',
                  'superkingdom': 'D',
                  'class': 'C',
                  'family': 'F',
                  'phylum': 'P'}

    print('... loading amps')
    df = pd.read_table('data/gmsc_meta_taxo.tsv.gz',
                       sep='\t',
                       header='infer',
                       low_memory=False)

    print('... loading reference')
    reftab = pd.read_table('data/conv_pgen_to_gtdb_list.tsv',
                           sep='\t',
                           header='infer',
                           low_memory=False)

    print('... processing pairs')
    newdf = []
    dfd = df[['level', 'name']]

    with tqdm(total=len(dfd),
              desc="Processed pairs",
              bar_format="{l_bar}{bar} [ time left: {remaining} ]") as pbar:
        for idx, l, taxon in dfd.itertuples():
            if l != 'no rank' and l == l:  # excludes no rand and n.a. taxonomy levels
                v = tab_levels[l]
                g = reftab[reftab[v] == taxon]['cluster'].drop_duplicates().tolist()
                if len(g) > 1 or len(g) == 0:
                    g = '*'               
            else:
                g = '*'
            newdf.append(g)
            pbar.update(1)
                
    print('... creating column')
    df['specI'] = newdf

    print('... saving table')
    df.to_csv('data/complete_amps_associated_taxonomy.tsv.gz',
              sep='\t', 
              header=True,
              index=None)     
              
