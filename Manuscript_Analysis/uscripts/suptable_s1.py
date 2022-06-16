import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

data = pd.read_table('data/samples-min500k-assembly-prodigal-stats.tsv.gz')
amps = pd.read_table('data/gmsc_amp_genes_envohr_source.tsv.gz')

namps = amps[['amp',
              'sample',
              'general_envo_name']]
              
namps = namps.drop_duplicates()
namps = namps.groupby('sample')
namps = namps.agg('size')
namps = namps.reset_index()

namps = namps.rename({0: 'amps',
                      'sample': 'sample_accession'},
                     axis=1)

a = data.merge(on='sample_accession',
               right=namps)
               
b = data[~data.sample_accession.isin(namps.sample_accession)]
data = pd.concat([a, b])
data.amps = data.amps.fillna(0)

data['amps_per_assembled_Mbp'] = data.amps * 1_000_000 / data.assembly_total_length

envo = pd.read_table('data/reduced_metadata.tsv.gz')
data = data.merge(on='sample_accession', right=envo)

## supp table S1
sup1 = data[['sample_accession', 'general_envo_name',
             'inserts_raw', 'assembly_total_length',
             'assembly_N50', 'prodigal_total_orfs',
             'smORFs', 'amps']].copy()

sup1.columns = ['sample', 'habitat', 'raw inserts',
                'assembled bp', 'N50', 'ORFs+smORFs',
                'smORFs', 'non-redundant AMPs']

sup1.to_csv('supplementary_table_S1.tsv', sep='\t', header=True, index=None)
