#!/usr/bin/env python

# Here we will show how c_AMPs from AMPSphere are relatively similar to those in the training set of Macrel (Santos-Júnior et al., 2020) and DRAMP 3.0 (Shi et al., 2021). To that, we will use pre-computed features calculated using Macrel's internal scripts from the peptide sequences in those three databases.

import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib import cm

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

def getlims(feat):
    feats = [
            ampsphere[feat].values,
            dramp[feat].values,
            macrel[feat].values
            ]
    xs = np.concatenate(feats)
    xs.sort()
    cut = 10
    return xs[cut], xs[-cut]

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

fig = plt.figure()
for ix,k in enumerate(graphkey):
    pos_i = ix % 3
    pos_j = ix // 3
    feat, label = graphkey[k]
    ax0 = fig.add_axes([.025 + pos_i * .3, 0.075 + pos_j * .3, .25, .12])
    ax1 = fig.add_axes([.025 + pos_i * .3, 0.1125 + pos_j * .3, .25, .12], sharex=ax0)
    ax2 = fig.add_axes([.025 + pos_i * .3, 0.150 + pos_j * .3, .25, .12], sharex=ax0)
    for data,data_label,ax,c in zip([ampsphere, dramp, macrel],
                                ['AMPSphere', 'DRAMP', 'Macrel'],
                                [ax0,ax1,ax2],
                                cm.Dark2.colors,
                                ):
        ax.patch.set_visible(False)
        sns.kdeplot(ax=ax,
                    fill=True,
                    bw_method='silverman',
                    cut=0,
                    bw_adjust=1.2,
                    data=data,
                    x=feat,
                    label=data_label,
                    color=c,
                    )
        ax.set_xlim(getlims(feat))
        ax.set_yticks([])
        ax.set_ylabel(None)
        ax.set_xlabel(None)
        ax.spines['left'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if k == 'a':
           ax.legend()
    for ax in [ax1, ax2]:
        ax.spines['bottom'].set_visible(False)
        ax.get_xaxis().set_visible(False)

    ax2.set_title(label,
                 fontfamily='Sans Serif',
                 fontsize='large',
                 loc='left')
