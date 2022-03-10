def analyze_level_III(analysis_folder):
    '''
    Loads fasta file for AMPSphere and 
    create a dataframe with the sequence,
    access and family; also computes
    the families which will be selected by
    size to next computations.
    '''
    import gzip
    import pandas as pd
    from collections import Counter
    from Bio import SeqIO
        
    # loading info about peptides and families
    infile = f'{analysis_folder}/AMPSphere_v.2021-03.faa.gz'
    lv3 = dict()
    for record in SeqIO.parse(gzip.open(infile, 
                                        'rt',
                                        encoding='utf-8'),
                              'fasta'):
        lv3[record.id] = [str(record.seq),
                          record.description.split(' | ')[1]]

    # selecting families of at least 8 peptides
    select_fam = [v[1] for k, v in lv3.items()]
    select_fam = [k for k, v in Counter(select_fam).items() if v >= 8]
    select_fam = sorted(select_fam)

    # converting dictionary to a dataframe
    lv3 = pd.DataFrame.from_dict(lv3,
                                 orient='index',
                                 columns=['sequence',
                                          'family'])
    lv3 = lv3.reset_index()
    lv3 = lv3.rename({'index': 'access'}, axis=1)

    return (lv3, select_fam)


def organize_house(analysis_folder):
    '''
    Create folders to receive the data generated
    by this script
    '''
    import os
    
    # creating folders for results
    os.makedirs(f'{analysis_folder}/families', exist_ok=True)
    subdirs = ['fastas',
               'aln',
               'tree_nwk',
               'tree_fig',
               'hmm',
               'hmm_logo']
    
    for folder in subdirs:
        os.makedirs(f'{analysis_folder}/families/{folder}',
                    exist_ok=True)


def fasta_files(lv3, select_fam, analysis_folder):
    '''
    Generates the fasta files per family
    '''
    import pandas as pd
    
    # computing results per family
    for f in select_fam:
        df = lv3[lv3['family'] == f]
        ofile = f'{analysis_folder}/families/fastas/{f}.faa'  
        with open(ofile, 'w') as db:
            for row in df.itertuples():
                db.write(f'>{row.access}\n{row.sequence}\n')


def alignments(select_fam, analysis_folder):
    '''
    Computes alignments from
    the fasta files generated
    before
    '''                
    from Bio.Align.Applications import MuscleCommandline
    
    for f in select_fam:
        ifile = f'{analysis_folder}/families/fastas/{f}.faa'
        ofile = f'{analysis_folder}/families/aln/{f}.aln'
        muscle_cline = MuscleCommandline(input=ifile,
                                         out=ofile,
                                         diags=True,
                                         maxiters=1)
        muscle_cline()
        

def trees(select_fam, analysis_folder):
    '''
    Computes trees and their ascII representation
    from the alignment files generated before
    '''                
    from .phylogen import treebuilder, draw_ascii

    for f in select_fam:
        # compute tree
        ofile = f'{analysis_folder}/families/tree_nwk/{f}.nwk'
        i2file = f'{analysis_folder}/families/aln/{f}.aln'
        treebuilder(i2file,
                    ofile,
                    'WAG+CAT',
                    1000)

        # compute image
        ifile = f'{analysis_folder}/families/tree_nwk/{f}.nwk'
        ofile = f'{analysis_folder}/families/tree_fig/{f}.ascii'
        draw_ascii(ifile,
                   ofile,
                   column_width=80)


def hmm(select_fam, data_folder, analysis_folder):
    '''
    Computes HMM profile and logo
    from the alignment files generated before
    '''                
    from .hmm import hbuild
    
    for f in select_fam:
        # computing HMM profile
        ofile = f'{analysis_folder}/families/hmm/{f}.hmm'
        i2file = f'{analysis_folder}/families/aln/{f}.aln'
        
        hbuild(i2file,
               ofile,
               f)


def hmmlogo(select_fam, data_folder, analysis_folder):
    '''
    Computes HMM logo for each family
    '''
    from .hmm import pict_hmmlogo
    
    for f in select_fam:
        # computing HMM logo
        alph = f'{data_folder}/alph.json'
        cmap = f'{data_folder}/cmap.json' 
        ifile = f'{analysis_folder}/families/hmm/{f}.hmm'
        ofile = f'{analysis_folder}/families/hmm_logo/{f}.svg'
        pict_hmmlogo(alph, cmap, ifile, ofile)


def process_cluster():
    import os
    
    data_folder = 'data/'
    analysis_folder = 'analysis/'
    
    for folder in [data_folder, analysis_folder]:
        os.makedirs(folder, exist_ok=True)
        
    print('Retrieve sequences and families')
    lv3, select_fam = analyze_level_III(analysis_folder)
    print('Generate folders')
    organize_house(analysis_folder)
    print('Generate fasta per family')
    fasta_files(lv3, select_fam, analysis_folder)
    print('Align peptides')
    alignments(select_fam, analysis_folder)
    print('Drawing trees')
    trees(select_fam, analysis_folder)
    print('Calculating Hidden Markov Models')
    hmm(select_fam, data_folder, analysis_folder)
    print('Calculating HMM logos')
    hmmlogo(select_fam, data_folder, analysis_folder)

