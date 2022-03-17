from utils.download_files import inputsgen
from utils.mundi_map import plot_mundi_map
from utils.homologs import homologs
from utils.henvo import heatmap_environments
from utils.hbody import heatmap_bodysites


def organize_files():
    import os
    flist = ['figure_1a_metagenomes_samples_distribution.svg',
             'panelB_homologs_search.txt',
             'figure_1c_heatmap_environments_overlap.svg',
             'figure_1d_heatmap_humanbodysites_overlap.svg']
    os.makedirs('analysis/', exist_ok=True)
    for f in flist:
        os.rename(f, f'analysis/{f}')


def main():
    print('Retrieving inputs')
    inputsgen()
    print('Generating figure 1A')
    plot_mundi_map()
    print('Generating data for figure 1B')
    homologs()
    print('Generating figure 1C')
    heatmap_environments()
    print('Generating figure 1D')
    heatmap_bodysites()

    
if __name__ == "__main__":
    main()
