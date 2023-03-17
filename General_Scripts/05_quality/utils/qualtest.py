def getref():
    '''
    Quick function to load file once
    '''
    import pandas as pd
    
    ref = pd.read_table('data/gmsc_amp_genes_envohr_source.tsv.gz',
                        sep='\t',
                        header='infer',
                        low_memory=False)
                        
    ref = ref[['gmsc', 'amp']]
    return ref
    
    
def filter_amps(supply, genes, testname, ref):
    '''
    Filter AMPs using a set of genes supposed to be
    either failing or passing the test (red or green)
    '''
    import pandas as pd
    
    all_amps = set(ref.amp)
    testresult = dict()

    if supply == 'red':
        red_amps = set(ref[ref.gmsc.isin(genes)].amp)
        green_amps = all_amps - red_amps        

    if supply == 'green':
        green_amps = set(ref[ref.gmsc.isin(genes)].amp)
        red_amps = all_amps - green_amps        

    for i in red_amps: testresult[i] = 'Failed'
    for i in green_amps: testresult[i] = 'Passed'

    df = pd.DataFrame.from_dict(testresult,
                                orient='index')
                                
    df = df.reset_index()
    df = df.rename({'index': 'AMP', 0: testname}, axis=1)
    
    return df    
    
    
def antifam_process(ref):
    '''
    Processes the results from ANTIFAM obtained for GMSC
    '''
    import gzip
    import pandas as pd
        
    red_genes = set(pd.read_table('data/antifam_100.tsv.gz',
                                  header=None)[0])

    df = filter_amps('red',
                     red_genes,
                     'Antifam',
                     ref)

    df.to_csv('antifam_list.tsv',
              sep='\t',
              header=True,
              index=None)
    

def terminal_place_process(ref):
    '''
    Processes the results from the 
    quality test of AMPs position in a contig
    spotting those AMPs happening in a terminal 
    position, which are possibily fragmented 
    genes
    '''
    import pandas as pd
    
    green_genes = pd.DataFrame()
    for chunk in pd.read_table('data/coordinates_test_passed.tsv.gz',
                               sep='\t',
                               header=None,
                               names=['gmsc', 'status'],
                               chunksize=10_000_000):
        chunk = chunk[chunk['status'] == 'T']
        chunk = chunk[chunk['gmsc'].isin(ref.gmsc)]
        green_genes = pd.concat([green_genes, chunk])

    del chunk
    green_genes = set(green_genes['gmsc'])
        
    df = filter_amps('green',
                     green_genes,
                     'Coordinates',
                     ref)

    df.to_csv('coordinates_list.tsv',
              sep='\t',
              header=True,
              index=None)


def metaP(ref):
    '''
    Processes the results from the 
    quality test of AMPs detection in metaproteomes
    '''
    import pandas as pd
    
    green_genes = set(pd.read_table('data/metaproteome_100.tsv.gz', header=None)[0])

    df = filter_amps('green',
                     green_genes,
                     'metaproteomes',
                     ref)

    df.to_csv('metaproteome_list.tsv',
              sep='\t',
              header=True,
              index=None)


def metaT(ref):
    '''
    Processes the results from the 
    quality test of AMPs detection via metatranscriptomics
    '''
    import pandas as pd

    green_amps = set()
    for chunk in pd.read_table('data/metaT_100.tsv.gz',
                               header=None,
                               chunksize=10_000_000):
        chunk = chunk[chunk[0].isin(ref.gmsc)][0]
        amps = set(ref[ref.gmsc.isin(chunk)].amp)
        green_amps = green_amps.union(amps)

    all_amps = set(ref.amp)
    red_amps = all_amps - green_amps

    testresult = dict()
    for i in red_amps: testresult[i] = 'Failed'
    for i in green_amps: testresult[i] = 'Passed'

    df = pd.DataFrame.from_dict(testresult,
                                orient='index')
    df = df.reset_index()
    df = df.rename({'index': 'AMP',
                    0: 'metatranscriptomes'},
                   axis=1)

    df.to_csv('metatranscriptome_list.tsv',
              sep='\t',
              header=True,
              index=None)


def assemble_tests():
    '''
    Assemble all results for separated quality tests
    '''
    import pandas as pd
    import matplotlib.pyplot as plt
    from collections import Counter

    print('... loading test results')
    d1 = pd.read_table('antifam_list.tsv', sep='\t', header='infer')
    d2 = pd.read_table('data/RNAcode_out_wlfam.tsv.xz', sep='\t', header='infer')
    d2.rename({'RNAcode_passed': 'RNAcode'}, axis=1, inplace=True)
    d3 = pd.read_table('metaproteome_list.tsv', sep='\t', header='infer')
    d4 = pd.read_table('metatranscriptome_list.tsv', sep='\t', header='infer')
    d5 = pd.read_table('coordinates_list.tsv', sep='\t', header='infer')
     
    print('... merging all tests')
    d6 = d1.merge(on='AMP', right=d2)
    d6 = d6.merge(on='AMP', right=d3)
    d6 = d6.merge(on='AMP', right=d4)
    d6 = d6.merge(on='AMP', right=d5)
    d6.sort_values(by='AMP', inplace=True)

    print('... returning experimental evidence of peptides')
    intdf = d6[['metaproteomes', 'metatranscriptomes']].replace('Passed', 1).replace('Failed', 0)
    intdf = intdf.sum(axis=1)
    intdf = intdf.replace([0, 1, 2],
                          ['Failed', 'Passed', 'Passed'])

    print('... counting amps per class')
    expevd = pd.DataFrame.from_dict(Counter(intdf), orient='index').T
    antifam = pd.DataFrame.from_dict(Counter(d6.Antifam), orient='index').T
    rna = pd.DataFrame.from_dict(Counter(d6.RNAcode), orient='index').T
    terminal = pd.DataFrame.from_dict(Counter(d6.Coordinates), orient='index').T

    print('... dataframing')
    df = pd.concat([antifam, terminal, expevd, rna])
    df.index = ['Antifam', 'Terminal placement', 'Experimental evidence', 'RNAcode']
    df = df.fillna(0)
    df = df[['Passed', 'Failed', 'Not tested']]

    print('... plot')
    df.plot.barh(stacked=True, cmap='Dark2')
    plt.xlabel('AMP candidates')
    plt.ylabel('Quality tests')
    plt.tight_layout()
    plt.savefig('figure_S1a_amp_quality.svg')

    print('... save files')
    d6.to_csv('quality_assessment.tsv', sep='\t', header=True, index=None)


def quality():
    print('Grabing reference')
    ref = getref()
    print('Processing Antifam results')
    antifam_process(ref)
    print('Processing coordinates results')
    terminal_place_process(ref)
    print('Processing metaproteomes results')
    metaP(ref)
    print('Processing metatranscriptomes results')
    metaT(ref)
    print('Assembling results')
    assemble_tests()
    
