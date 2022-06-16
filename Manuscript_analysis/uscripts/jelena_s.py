import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

print('# load datasets')
ampsphere = pd.read_csv('data/ampsphere_v2022-03.features.tsv.gz', sep="\t")
dramp = pd.read_csv('data/dramp_v3.features.tsv.gz', sep="\t")
macrel = pd.read_csv('data/macrel_trainpos.features.tsv.gz', sep="\t")


print('# calculating length')
ampsphere['length'] = ampsphere.sequence.apply(lambda x: len(x))
dramp['length'] = dramp.sequence.apply(lambda x: len(x))
macrel['length'] = macrel.sequence.apply(lambda x: len(x))

ampsphere.drop(['group',
                'sequence'],
               axis=1,
               inplace=True)

dramp.drop(['group',
            'sequence'],
           axis=1,
           inplace=True)

macrel.drop(['group',
             'sequence'],
            axis=1,
            inplace=True)

print('# plotting and subplotting (3 rows x 3 columns for 9 features)')
fig, axarr = plt.subplot_mosaic([
                                 ['a)', 'b)', 'c)'],
                                 ['d)', 'e)', 'f)'],
                                 ['g)', 'h)', 'i)'],
                                 ])

def getminmax(feat):
    a, b = [], []
    for j in [ampsphere, dramp, macrel]:
        a.append(j[feat].min())
        b.append(j[feat].max())
    return min(a), max(b)

print('# setting keys')
graphkey = {
            'a)': ['length', 'Length (residues)'], 
            'b)': ['smallAA', 'Small residues'],
            'c)': ['basicAA', 'Basic residues'],
            'd)': ['pI', 'Isoelectric point'],
            'e)': ['charge', 'Charge at pH 7.0'],
            'f)': ['aindex', 'Aliphatic index'],
            'g)': ['instaindex', 'Instability index'],
            'h)': ['boman', 'Boman index'],
            'i)': ['hmoment', 'Hydrophobic moment']
            }
            
for k in graphkey:
   feat, label = graphkey[k]
   print(f'# plotting pannel {k} - {feat}')
   sns.kdeplot(ax=axarr[k],
               data=ampsphere,
               x=feat,
               label='AMPSphere')
   sns.kdeplot(ax=axarr[k],
               data=dramp,
               x=feat,
               label='DRAMP')
   sns.kdeplot(ax=axarr[k],
               data=macrel,
               x=feat,
               label='Macrel')
   axarr[k].set_xlabel(label)
   axarr[k].set_xlim(getminmax(feat))
   if k == 'a)':
       axarr[k].legend()
   if k not in ['a)', 'd)', 'g)']:
       axarr[k].set_yticks([])
       axarr[k].set_ylabel(None)

for label, ax in axarr.items():
    ax.set_title(label,
                 fontfamily='Sans Serif',
                 fontsize='large',
                 loc='left')

plt.tight_layout()
plt.savefig('all_features.svg')

