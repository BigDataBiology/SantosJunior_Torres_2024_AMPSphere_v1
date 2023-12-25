#!/usr/bin/env python

# c_AMPs were converted to a 8-letter alphabet, then were hierarchically clustered using CD-Hit at three levels: 100%, 85% and 75% of identity. To make the validation, we sampled 1000 sequences in triplicate, aligned them against the cluster representative and calculated the E-value of the alignment.

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

lv1 = pd.read_table('data/output_clustering_significance_levelI.tsv.gz')
lv2 = pd.read_table('data/output_clustering_significance_levelII.tsv.gz')
lv3 = pd.read_table('data/output_clustering_significance_levelIII.tsv.gz')


# preparing variables
lv1['min_id'] = lv1[['identity', 'gap_identity']].min(axis=1)
lv2['min_id'] = lv2[['identity', 'gap_identity']].min(axis=1)
lv3['min_id'] = lv3[['identity', 'gap_identity']].min(axis=1)

lv1['Log(E-value)'] = np.log10(lv1.evalue)
lv2['Log(E-value)'] = np.log10(lv2.evalue)
lv3['Log(E-value)'] = np.log10(lv3.evalue)


# the proportion of significant alignments at each cluster level
pct1 = lv1.eval('`Log(E-value)` <= -5').sum() / len(lv1)
pct2 = lv2.eval('`Log(E-value)` <= -5').sum() / len(lv2)
pct3 = lv3.eval('`Log(E-value)` <= -5').sum() / len(lv3)

level_pct = [
    ('I', pct1),
    ('II', pct2),
    ('III',pct3),
    ]


# plotting
fig, axarr = plt.subplot_mosaic([['a)', 'b)', 'c)']],
                                constrained_layout=True)

sns.scatterplot(ax=axarr['a)'],
                data=lv1,
                x='min_id',
                y='Log(E-value)',
                hue='replicate',
                s=2.5,
                alpha=0.5,
                palette='Dark2')

sns.scatterplot(ax=axarr['b)'],
                data=lv2,
                x='min_id',
                y='Log(E-value)',
                hue='replicate',
                s=2.5, alpha=0.5,
                legend=False,
                palette='Dark2')

sns.scatterplot(ax=axarr['c)'],
                data=lv3,
                x='min_id',
                y='Log(E-value)',
                hue='replicate',
                s=2.5,
                alpha=0.5,
                legend=False,
                palette='Dark2')

for idx, (label, ax) in enumerate(axarr.items()):
    level, pct = level_pct[idx]
    ax.set_title(f'Clstr. Lv. {level} [{pct:.1%}]', fontfamily='Sans Serif', fontstyle='italic')
    ax.axhline(y=-5., color='black', linestyle='dashed', linewidth=1.0)
    ax.set_xlabel('Identity (%)')
    ax.set_xlim(0, 100)
    ax.set_ylim(20, -80)
    if idx > 0:
        ax.set_ylabel(None)
        ax.set_yticks([])
    if idx != 1:
        ax.set_xlabel(None)

fig.savefig('data/output_clustering_validation.svg')

