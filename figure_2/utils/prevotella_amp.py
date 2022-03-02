def prevotella():
    '''
    Generates a table with the Prevotella species and
    the number of AMPs for each, also includes the
    host and environment information -- it will be used
    as metadata for the tree drawing in the iToL web
    application
    '''
    import pandas as pd
    import matplotlib.pyplot as plt

    # the file "taxonomy_annotation.tsv" is generated during code for Fig. 2a
    data = pd.read_table('taxonomy_annotation.tsv')
    
    # clean data and select only prevotella
    # species already isolated and described
    ds = data[['amp', 'fixed', 'name']]
    ds = ds[(ds.name.str.contains('Prevotella ')) & (~ds.name.str.contains('Prevotella sp'))]
    ds = ds.sort_values(by='name')

    # prevotella_species_list was manually curated and is available after input download
    pspecies = pd.read_table('data/prevotella_species_list.tsv')
    pspecies_list = pspecies['Species name'].str.replace('P. ', 'Prevotella ').tolist()

    # counts prevotella unique AMPs
    pspecies_amps = []
    for ps in pspecies_list:
        amps = len(ds[ds.name.str.contains(ps)]['amp'].drop_duplicates())
        print(f'{ps}\t{amps}')
        pspecies_amps.append(amps)

    # add prevotella AMP counts to the dataframe
    pspecies['AMPs'] = pspecies_amps
    
    # saving text table to perform the tree in the iTOL website
    pspecies.to_csv('prevotella_species_amp_counts.tsv',
                    sep='\t',
                    header=True,
                    index=None)
                    

def alignment():
    """
    Perform the alignment of 16S sequences from Prevotella genus.
    """
    from Bio.Align.Applications import MafftCommandline
    
    mafft_cline = MafftCommandline(input='data/37185.fasta')
    print(mafft_cline)

    stdout, stderr = mafft_cline()

    with open ('37185.aln','w') as handle:
        handle.write(stdout)


def tree():
    '''
    Calculate the Newick tree for Prevotella genus
    starting from the alignment of the 16S gene
    '''
    import subprocess
    
    alignment()
    subprocess.call(['fasttree', 
                     '-nt',
                     '-gtr',
                     '-pseudo',
                     '1000',
                     '-out',
                     '37185.nwk',
                     '37185.aln'])


def prevamp():
    print('Computes the matrix of AMPs per Prevotella species')
    prevotella()
    print('Generate tree for Prevotella genus')
    tree()
    
