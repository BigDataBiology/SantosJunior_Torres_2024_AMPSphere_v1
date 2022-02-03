import pandas as pd
import matplotlib.pyplot as plt

# first execute codes for Fig.1a
# the file "taxonomy_annotation.tsv" will be generated and used in this script
data = pd.read_table('../figure_2a/taxonomy_annotation.tsv')
ds = data[['amp', 'fixed', 'name']]
ds = ds[(ds.name.str.contains('Prevotella ')) & (~ds.name.str.contains('Prevotella sp'))]
ds = ds.sort_values(by='name')

# table available in > ubuntu@aws.big-data-biology.org:/share/work/Celio/files_for_figures/genes_origins/
pspecies = pd.read_table('prevotella_species_list.tsv')
pspecies_list = pspecies['Species name'].str.replace('P. ', 'Prevotella ').tolist()

pspecies_amps = []
for ps in pspecies_list:
    amps = len(ds[ds.name.str.contains(ps)]['amp'].drop_duplicates())
    print(f'{ps}\t{amps}')
    pspecies_amps.append(amps)

pspecies['AMPs'] = pspecies_amps
pspecies.to_csv('prevotella_species_amp_counts.tsv',
                sep='\t',
                header=True,
                index=None)
