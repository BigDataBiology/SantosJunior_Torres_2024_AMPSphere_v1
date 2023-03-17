from utils.amp_density_gut_traveling import density_travel
#from utils.download_files import inputsgen
from utils.taxonomy import taxon_analysis
from utils.prevotella_amp import prevamp
from utils.family_sizes import family_size
from utils.core import plot_core
from utils.itol_gtdb import itol_gtdb_prep

def organize_files():
    '''
    Realocate outputs into analysis folder
    '''
    import os
    import glob
    
    flist = ['37185.aln',
             'result_tax.txt',
             '37185.nwk',
             'taxonomy_annotation.tsv',
             'amps_all.count_core.tsv',
             'families_all.count_core.tsv',
             'prevotella_species_amp_counts.tsv']

    os.makedirs('analysis/figures/', exist_ok=True)
    for f in flist:
        os.rename(f, f'analysis/{f}')
    for f in glob.glob('*.svg'):
        os.rename(f, f'analysis/figures/{f}')
    
    
def main():
#    print('Downloading inputs')
#    inputsgen()
    print('Generating figure 2A')
    taxon_analysis()
    print('Generating figure 2B')
    # this generates the two files needed
    # as inputs for the iToL web application
    # just enter the tree and the table as 
    # annotation and the tree from the 
    # pannel B is regenerated
    prevamp()
    print('Generating figure 2C')
    family_size()
    print('Generating figure 2D')
    plot_core()
    print('Organizing results')
    organize_files()
    print('Testing density in gut travelling microbes')
    density_travel()
    print('Preparing tree with densities')
    itol_gtdb_prep()
    
    
if __name__ == "__main__":
    main()
