from utils.qualtest import quality
from utils.ugenes import ugenes_plot
from utils.AMPfeatures import amplen
#from utils.download_files import inputsgen
from utils.enrich import enrichment_analysis

def organize_outputs():
    '''
    Generate analysis/ folder and organize the outputted results
    into it
    '''
    import os
    
    os.makedirs('analysis/', exist_ok=True)
    os.makedirs('analysis/figures/', exist_ok=True)

    flist = ['antifam_list.tsv',
             'coordinates_list.tsv',
             'metaproteome_list.tsv',
             'metatranscriptome_list.tsv',
             'quality_assessment.tsv',
             'quality_candidates.txt',
             'quality_families.txt',
             'rnacode_list.tsv',
             'high_quality_candidates.txt']

    figlist = ['figure_S1a_amp_quality.svg',
               'figure_S1b_unique_genes.svg',
               'figure_S1c_histogram_AMPs_length.svg']

    for fin in flist: os.rename(fin, f'analysis/{fin}')
    for fin in figlist: os.rename(fin, f'analysis/figures/{fin}')
    
    
def main():
#    print('Downloading inputs')
#    inputsgen()
    print('Generating figure S1A')
    quality()
    print('Testing statistics')
    enrichment_analysis()
    print('Generating figure S1B')
    ugenes_plot()
    print('Generating figure S1C')
    amplen()
    print('Generating figure S1D')
    heatmap_mammal_guts()
    print('Organize files')
    organize_outputs()
    
    
if __name__ == "__main__":
    main()
