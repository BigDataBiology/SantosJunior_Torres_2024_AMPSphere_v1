from utils.mundi_map import plot_mundi_map
from utils.homologs import homologs
from utils.henvo import heatmap_environments
from utils.hbody import heatmap_bodysites
from utils.download_files import inputsgen

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
