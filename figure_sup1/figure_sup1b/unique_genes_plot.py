import lzma
import pandas as pd
import numpy as np

from Bio import SeqIO
from collections import Counter


# File 'AMPSphere_v.2021-03.fna.xz' available in AMPSphere Zenodo repository

# create a table with sequence and AMP headers
headers, seqs = [], []
with lzma.open('AMPSphere_v.2021-03.fna.xz', 'rt') as infile:
    for record in SeqIO.parse(infile, 'fasta'):
        headers.append(record.description.split()[1])
        seqs.append(str(record.seq))

df = pd.DataFrame(np.array([headers, seqs]).T, columns=['AMP', 'gene'])
df.AMP = [x[1] for x in df.AMP]
df = df.drop_duplicates()

ds = pd.DataFrame.from_dict(Counter(dict(Counter(df.AMP)).values()),
                            orient='index',
                            columns=['number_of_amps'])
ds = ds.reset_index()
ds = ds.rename({'index': 'number_of_ugenes'},
               axis=1)
               
ds = ds.sort_values(by='number_of_ugenes')
amps = ds[ds.number_of_ugenes >= 3]['number_of_amps'].sum()

ds = ds.set_index('number_of_ugenes')
ds = pd.DataFrame([ds.loc[1]['number_of_amps'],
                   ds.loc[2]['number_of_amps'],
                   amps],
                  index=['1', '2', '>3'])

ds.plot.bar(legend=False, cmap='Dark2')
plt.xlabel('Number of unique genes')
plt.ylabel('AMP candidates')
plt.savefig('unique_genes_distribution.svg')

