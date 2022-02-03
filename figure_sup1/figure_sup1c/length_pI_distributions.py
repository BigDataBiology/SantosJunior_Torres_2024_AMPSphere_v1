import pandas as pd
import matplotlib.pyplot as plt

# loading features
# file available in > ubuntu@aws.big-data-biology.org:/share/work/Celio/files_for_figures/features/
data = pd.read_table('AMPsphere_features.tsv', sep='\t', header='infer')
data['length'] = data.sequence.str.len()

# preparing the histogram of isolectric points
data.length.hist(bins=100, grid=False)
plt.xlim(5,100)
plt.ylabel('AMP candidates')
plt.xlabel('Length in residues')
plt.savefig('amp_candidates_length.svg')
plt.close()

# preparing the histogram of isolectric points
data.pI.hist(bins=100, grid=False)
plt.xlim(0, 15)
plt.ylabel('AMP candidates')
plt.xlabel('Isoelectric point')
plt.savefig('amp_candidates_pI.svg')
plt.close()

