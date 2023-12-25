#!/usr/bin/env python

# Here we will show how c_AMPs from AMPSphere are relatively similar to those in the training set of Macrel (Santos-Júnior et al., 2020) and DRAMP 3.0 (Shi et al., 2021). To that, we will use pre-computed features calculated using Macrel's internal scripts from the peptide sequences in those three databases.

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

ampsphere = pd.read_csv('../data_folder/ampsphere_v2022-03.features.tsv.gz', sep="\t")
dramp = pd.read_csv('../data_folder/dramp_v3.features.tsv.gz', sep="\t")
macrel = pd.read_csv('../data_folder/macrel_trainpos.features.tsv.gz', sep="\t")

# Macrel's internal script to compute features did not use the protein length
# as a feature, and therefore, to compare it through these different datasets,
# we will need to compute it.
# We also remove columns that Macrel uses in its internal pipeline and we won't
# be using here.

for feats in [ampsphere, dramp, macrel]:
    feats['length'] = feats.sequence.str.len()
    feats.drop(['group',
            'sequence'],
           axis=1,
           inplace=True)

def getminmax(feat, ampsphere, dramp, macrel):
    a, b = [], []
    for j in [ampsphere, dramp, macrel]:
        a.append(j[feat].min())
        b.append(j[feat].max())
    return min(a), max(b)


graphkey = {
            'a': ['length',     'Length (residues)'],
            'b': ['smallAA',    'Small residues'],
            'c': ['basicAA',    'Basic residues'],
            'd': ['pI',         'Isoelectric point'],
            'e': ['charge',     'Charge at pH 7.0'],
            'f': ['aindex',     'Aliphatic index'],
            'g': ['instaindex', 'Instability index'],
            'h': ['boman',      'Boman index'],
            'i': ['hmoment',    'Hydrophobic moment']
            }

fig, axarr = plt.subplot_mosaic([
                                 ['a', 'b', 'c'],
                                 ['d', 'e', 'f'],
                                 ['g', 'h', 'i'],
                                 ])


for k in graphkey:
    feat, label = graphkey[k]
    ax = axarr[k]
    ax.clear()
    for data,data_label in zip([ampsphere, dramp, macrel],
                               ['AMPSphere', 'DRAMP', 'Macrel']):
        sns.kdeplot(ax=ax,
                    fill=True,
                    bw_method='silverman',
                    cut=0,
                    bw_adjust=1.2,
                    data=data,
                    x=feat,
                    label=data_label,
                    )
    ax.set_xlim(getminmax(feat,
                            ampsphere,
                            dramp,
                            macrel))
    if k == 'a':
       ax.legend()

    ax.set_yticks([])
    ax.set_ylabel(None)
    ax.set_title(label,
                 fontfamily='Sans Serif',
                 fontsize='large',
                 loc='left')

fig.tight_layout()
