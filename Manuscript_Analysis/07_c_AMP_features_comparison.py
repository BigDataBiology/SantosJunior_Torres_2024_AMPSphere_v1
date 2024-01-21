# Here we will show how c_AMPs from AMPSphere are relatively similar to those in the training set of Macrel (Santos-Júnior et al., 2020) and DRAMP 3.0 (Shi et al., 2021). To that, we will use pre-computed features calculated using Macrel's internal scripts from the peptide sequences in those three databases.

import pandas as pd
import numpy as np
import seaborn as sns
from scipy import stats
from matplotlib import pyplot as plt
from matplotlib import cm
from modlamp.descriptors import GlobalDescriptor
from modlamp.descriptors import PeptideDescriptor
from statsmodels.stats.multitest import multipletests

from os import makedirs
makedirs('figures', exist_ok=True)
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['xtick.labelsize'] = 7

ampsphere = pd.read_csv('../data_folder/ampsphere_v2022-03.features.tsv.gz', sep="\t")
dramp = pd.read_csv('../data_folder/dramp_v3.features.tsv.gz', sep="\t")
macrel = pd.read_csv('../data_folder/macrel_trainpos.features.tsv.gz', sep="\t")
neg_macrel = pd.read_csv('new_data/macrel_trainneg.features.tsv.gz', sep="\t")

# Macrel's internal script to compute features did not use the protein length
# as a feature, and therefore, to compare it through these different datasets,
# we will need to compute it.
# We also remove columns that Macrel uses in its internal pipeline and we won't
# be using here.

for feats in [ampsphere, dramp, macrel, neg_macrel]:
    feats['length'] = feats.sequence.str.len()
    gc = GlobalDescriptor(feats.sequence.tolist())
    gc.hydrophobic_ratio()
    feats['hydro_ratio'] = gc.descriptor.ravel()

    pepd = PeptideDescriptor(feats.sequence.tolist(), 'flexibility')
    pepd.calculate_global()
    feats['flexibility'] = pepd.descriptor.ravel()

    pepd = PeptideDescriptor(feats.sequence.tolist(), 'gravy')
    pepd.calculate_global()
    feats['gravy'] = pepd.descriptor.ravel()

    feats.drop(['group',
            'sequence'],
           axis=1,
           inplace=True)

def getlims(feat):
    xmin = ampsphere[feat].max()
    xmax = ampsphere[feat].min()
    for feats in [
            ampsphere[feat].values,
            dramp[feat].values,
            macrel[feat].values,
            neg_macrel[feat].values,
            ]:
        xs = feats.copy()
        xs.sort()
        cut = len(xs) // 100
        if xs[cut] < xmin:
            xmin = xs[cut]
        if xs[-cut] > xmax:
            xmax = xs[-cut]
    return xmin, xmax

for f in ['smallAA', 'basicAA', 'polarAA', 'aromaticAA', 'chargedAA']:
    for feats in [ampsphere, dramp, macrel, neg_macrel]:
        feats[f] *= 100
panels = [
          ('length',     'Length (residues)'),
          ('smallAA',    'Small residues (%)'),
          ('basicAA',    'Basic residues (%)'),
          ('polarAA',    'Polar residues (%)'),
          ('aromaticAA', 'Aromatic residues (%)'),
          ('chargedAA',  'Charged residues (%)'),
          ('pI',         'Isoelectric point'),
          ('charge',     'Charge at pH 7.0'),
          ('aindex',     'Aliphatic index'),
          ('instaindex', 'Instability index'),
          ('boman',      'Boman index'),
          ('hmoment',    'Hydrophobic\nmoment'),
          ('hydro_ratio','Hydrophobic ratio'),
          ('hydrophobicity','Hydrophobicity'),
          ('flexibility','Flexibility'),
          ('gravy',      'GRAVY'),
          ]

PANEL_WIDTH = 0.2
SUBPANEL_HEIGHT = 0.05

fig = plt.figure()
for ix,(feat,label) in enumerate(panels):
    pos_i = ix % 4
    pos_j = 3 - ix // 4
    ax0 = fig.add_axes([.025 + pos_i * .25, 0.075 + pos_j * .23, PANEL_WIDTH, SUBPANEL_HEIGHT])
    ax1 = fig.add_axes([.025 + pos_i * .25, 0.100 + pos_j * .23, PANEL_WIDTH, SUBPANEL_HEIGHT], sharex=ax0)
    ax2 = fig.add_axes([.025 + pos_i * .25, 0.125 + pos_j * .23, PANEL_WIDTH, SUBPANEL_HEIGHT], sharex=ax0)
    ax3 = fig.add_axes([.025 + pos_i * .25, 0.150 + pos_j * .23, PANEL_WIDTH, SUBPANEL_HEIGHT], sharex=ax0)
    for data,data_label,ax,c in zip([neg_macrel, macrel, dramp, ampsphere],
                                ['Macrel (neg)', 'Macrel (pos)', 'DRAMP', 'AMPSphere'],
                                [ax0, ax1, ax2, ax3],
                                cm.Dark2.colors[:4][::-1], # Use the first 4 colours, but have the negative group by red(ish)
                                ):
        ax.patch.set_visible(False)
        sns.kdeplot(ax=ax,
                    fill=True,
                    bw_method='silverman',
                    bw_adjust=(2. if feat.endswith('AA') else 1.),
                    cut=0,
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
        if ix == 0:
           ax.legend()
    for ax in [ax1, ax2, ax3]:
        ax.spines['bottom'].set_visible(False)
        ax.get_xaxis().set_visible(False)

    ax3.set_title(label,
                 fontfamily='Sans Serif',
                 fontsize='small',
                 loc='left')

fig.savefig('figures/07_c_AMP_features_comparison.svg', bbox_inches='tight')
fig.savefig('figures/07_c_AMP_features_comparison.png', bbox_inches='tight', dpi=300)

data_groups = [neg_macrel, macrel, dramp, ampsphere]
data_names = ['Macrel (neg)', 'Macrel (pos)', 'DRAMP', 'AMPSphere']
comparisons = []
for i in range(len(data_groups)):
    for j in range(i+1, len(data_groups)):
        for feat,label in panels:
            u,p = stats.mannwhitneyu(data_groups[i][feat], data_groups[j][feat])
            comparisons.append((data_names[i], data_names[j], feat, p))
            print(f'{data_names[i]} vs {data_names[j]}: {label}: p={p:.2e}')
comparisons = pd.DataFrame(comparisons, columns=['Group 1', 'Group 2', 'Feature', 'P-value'])
_, comparisons['padj'], _, _ = multipletests(comparisons['P-value'])
comparisons.to_csv('outputs/07_c_AMP_features_comparison.tsv', sep='\t', header=True, index=None)

