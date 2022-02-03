import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt

def classificationprop(x):
    if x >= 90: return 'core'
    if x >= 50: return 'shell'
    if x < 50: return 'accessory'

# File proGenomes2.1_specI_clustering.tab can be downloaded as:
# wget https://progenomes.embl.de/data/proGenomes2.1_specI_clustering.tab
df = pd.read_table('proGenomes2.1_specI_clustering.tab', sep='\t', header=None)
spec_size = dict(Counter(df[0]))

# File SPHERE_v.2021-03.levels_assessment.tsv.gz is available in the AMPSphere Zenodo repository
spheres = pd.read_table('SPHERE_v.2021-03.levels_assessment.tsv.gz', sep='\t', header='infer')
spheres = spheres[['AMP accession', 'SPHERE_fam level III']]
DD = spheres.set_index('AMP accession').to_dict()['SPHERE_fam level III']

# file complete_gmsc_pgenomes_metag.tsv is available in ubuntu@aws.big-data-biology.org:/share/work/Celio/files_for_figures/genes_origins/
refdata = pd.read_table('complete_gmsc_pgenomes_metag.tsv', sep='\t', header='infer')
refdata = refdata[refdata.is_metagenomic == False]
refdata = refdata[refdata.specI != '*']
refdata = refdata[['amp','sample','specI']]
refdata = refdata.sort_values(by=['amp', 'sample'])
refdata = refdata.drop_duplicates()
refdata['families'] = [DD[x] for x in refdata.amp]

## analyzing families
fams = refdata[['families', 'sample', 'specI']]
fams = fams.drop_duplicates()
fams = fams.drop('sample', axis=1)
fams = fams.sort_values(by=['families', 'specI'])

# counting clusters per family
ofile = open('families_all.count_core.tsv', 'w')
ofile.write('family\tspecI\tcounts\ttotal\tproportion\tclassification\n')

for df in fams.groupby('families'):
    f = df[0]
    D = dict(Counter(df[1].specI))
    for k, v in D.items():
        if spec_size[k] >= 10:
            nv = v * 100 / spec_size[k]
            total = spec_size[k]
            ofile.write(f'{f}\t{k}\t{v}\t{total}\t{nv}\t{classificationprop(nv)}\n')

ofile.close()

## analyzing amps
amps = refdata[['amp', 'sample', 'specI']]
amps = amps.drop_duplicates()
amps = amps.drop('sample', axis=1)

# counting clusters per family
ofile = open('amps_all.count_core.tsv', 'w')
ofile.write('amp\tspecI\tcounts\ttotal\tproportion\tclassification\n')

for df in amps.groupby('amp'):
    f = df[0]
    D = dict(Counter(df[1].specI))
    for k, v in D.items():
        if spec_size[k] >= 10:
            nv = v * 100 / spec_size[k]
            total = spec_size[k]
            ofile.write(f'{f}\t{k}\t{v}\t{total}\t{nv}\t{classificationprop(nv)}\n')

ofile.close()

