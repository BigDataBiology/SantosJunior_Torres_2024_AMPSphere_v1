from utils.qualtest import quality
from utils.ugenes import ugenes_plot
from utils.AMPfeatures import amplen
from utils.hmammal_guts import heatmap_mammal_guts
from utils.download_files import inputsgen
from utils.enrich import enrichment_analysis

def main():
    print('Downloading inputs')
    inputsgen()
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
    
if __name__ == "__main__":
    main()
