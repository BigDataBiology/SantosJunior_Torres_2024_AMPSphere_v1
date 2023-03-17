import pandas as pd
from glob import glob

fams, seqs = [], []
for f in glob('analysis/*_proteins.tsv.gz'):
    cluster = f.replace('_proteins.tsv.gz', '')
    cluster = cluster.replace('analysis/', '')
    df = pd.read_table(f, sep='\t', header='infer')
    n50 = df.loc[df.cutoff_prevalence == 50, 'percent of unique proteins'].tolist()[0]
    n95 = df.loc[df.cutoff_prevalence == 95, 'percent of unique proteins'].tolist()[0]
    ac = 100-n50
    shell = n50-n95
    seqs.append((cluster, ac, shell, n95))

seqs = pd.DataFrame(seqs,
                    columns=['cluster', 'accessory', 'shell', 'core'])

seqs.to_csv('progenomes_Core_Shell_Acc.tsv.gz',
            sep='\t',
            header=True,
            index=None)

for f in glob('analysis/*_families.tsv.gz'):
    cluster = f.replace('_families.tsv.gz', '')
    cluster = cluster.replace('analysis/', '')
    df = pd.read_table(f, sep='\t', header='infer')
    n50 = df.loc[df.cutoff_prevalence == 50, 'percent of families'].tolist()[0]
    n95 = df.loc[df.cutoff_prevalence == 95, 'percent of families'].tolist()[0]
    ac = 100-n50
    shell = n50-n95
    fams.append((cluster, ac, shell, n95))

fams = pd.DataFrame(fams,
                    columns=['cluster', 'accessory', 'shell', 'core'])

fams.to_csv('progenomes_fam_Core_Shell_Acc.tsv.gz',
            sep='\t',
            header=True,
            index=None)

