from utils.download_files import inputsgen
from utils.taxonomy import taxon_analysis
from utils.prevotella_amp import prevamp
from utils.family_sizes import family_size
from utils.core import plot_core

def main():
    print('Downloading inputs')
    inputsgen()
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
    
if __name__ == "__main__":
    main()
