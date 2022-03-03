from utils.download_files import inputsgen
from utils.gtdb2pgenomes import conversion_from_GTDB_to_pGenomes
from utils.preprocess import preprocess
from utils.progenomes_genes import ampsphere2progenomes
from utils.metadata import metadata
from utils.addspecI import addspecI
from utils.hr_envo import hrenvo

def main():
    print('Download inputs')
    inputsgen()
    print('Create conversion tables from GTDB into ProGenomes')
    conversion_from_GTDB_to_pGenomes()
    print('Linking taxonomy and GMSC genes')
    preprocess()
    print('Generating table of genes from ProGenomes2')
    ampsphere2progenomes()
    print('Processing metadata')
    metadata()
    print('Adding specI cluster info to the tables')
    addspecI()
    print('Adding metadata info and consolidating data')
    hrenvo()
    
