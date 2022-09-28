def anno(fin, ofile, dramp, bout):
    '''
    Add description to the dramp targets and save the final file
    '''
    import pandas as pd
    from .dramp_anno import dramp_anno
    from .utils import call_mmseqs
     
    dramp_anno = pd.DataFrame.from_dict(dramp_anno, orient='index')
    dramp_anno = dramp_anno.reset_index()
    dramp_anno = dramp_anno.rename({'index': 'target',
                                    0: 'description',
                                    1: 'activity',
                                    2: 'phase_of_testing',
                                    3: 'reference'},
                                   axis=1)
    
    print ('Searching reminiscent singletons against DRAMP, it can take a while')
    Nthreads = 3
    call_mmseqs(fin, bout, dramp, Nthreads, tmp='tmp', stdout=None)
        
    data = pd.read_table(bout,
                         sep='\t',
                         header=None,
                         names=['query', 'target', 'fident',
                                'alnlen', 'mismatch', 'gapopen',
                                'qstart', 'qend', 'tstart', 'tend',
                                'evalue', 'bits'])
 
    parsed = data[(data.fident >= 0.75) & (data.evalue <= 1e-5)]
    
    parsed = parsed.sort_values(by=['bits', 'evalue', 'fident'],
                                ascending=[False, True, False])
        
    parsed = parsed.merge(on='target', right=dramp_anno)
    
    parsed.to_csv(ofile,
                  sep='\t', header=True, index=None)
 

def dramp_anno():
    '''
    Annotates AMPSphere with DRAMP database
    '''
    import os
    
    data_folder = 'data/'
    analysis_folder = 'analysis/'
    
    for i in [data_folder, analysis_folder]:
        os.makedirs(i, exist_ok=True)

    # database to be searched against
    dramp = f'{data_folder}/DRAMP_GCP2019.fasta'
    # fasta input        
    fin = f'{analysis_folder}/AMPSphere_v.2022-03.faa.gz'
    # intermediate mmseqs result
    bout = f'{analysis_folder}/DRAMP_annotation.raw.tsv'
    # mmseqs parsed table outputted
    ofile = f'{analysis_folder}/DRAMP_anno_AMPSphere_v.2022-03.parsed.tsv.gz'

    anno(fin, ofile, dramp, bout)
    
