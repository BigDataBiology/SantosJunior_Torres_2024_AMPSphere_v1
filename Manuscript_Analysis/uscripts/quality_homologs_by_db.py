import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def getproportions(amps, quality):
    '''
    Get percent of passing AMPs per test
    in a given set of AMPs
    
    :inputs:
    amps - list of amps
    quality - pandas dataframe with the
              quality assessment per AMP
    
    :output:
    dataframe with the index as different quality
    tests and the percent of homologs passing the
    tests
    '''
    dftmp = quality[quality.AMP.isin(amps)]
    dftmp = dftmp.drop('AMP', axis=1)
    dftmp = dftmp.replace('Passed', 1)
    dftmp = dftmp.replace('Failed', 0)
    dftmp = dftmp.replace('Not tested', 0)
    return dftmp.sum(axis=0) * 100 / len(amps)
    
    
# load info
import pandas as pd

dramp = pd.read_table("data/dramp_candidates.txt.gz", header=None)[0].tolist()
gmgc = pd.read_table("data/gmgc_candidates.txt.gz", header=None)[0].tolist()
smprot = pd.read_table("data/SmProt_candidates.txt.gz", header=None)[0].tolist()
starpep = pd.read_table("data/starPepDB_candidates.txt.gz", header=None)[0].tolist()
sts = pd.read_table("data/STsORFs_candidates.txt.gz", header=None)[0].tolist()
qual = pd.read_table('data/quality_assessment.tsv.gz')

dramp = getproportions(dramp, qual)
gmgc = getproportions(gmgc, qual)
smprot = getproportions(smprot, qual)
starpep = getproportions(starpep, qual)
sts = getproportions(sts, qual)

df = pd.concat([dramp, gmgc, smprot, starpep, sts], axis=1)
df.columns = ['DRAMP v3.0', 'GMGC v1', 'SmProt 2', 'starPep45k', 'STsORFs']
df.to_csv('quality_by_db.tsv', sep='\t', header=True, index=True)
df = df.reset_index().rename({'index': 'Quality test'}, axis=1)
ampsphere_limits = getproportions(qual.AMP.tolist(), qual).to_dict()

x = df.melt(id_vars=['Quality test'], value_vars=df.columns[1:])

x.rename({'value': '% of homologs',
          'variable': 'Database'},
         axis=1,
         inplace=True)

# plot
fig, ax = plt.subplots()
sns.barplot(x='Quality test', y='% of homologs', hue='Database', data=x)

# iterate over range of number of rows
for i, c in enumerate(ampsphere_limits):
    ax.hlines(y=ampsphere_limits[c], xmin=i-0.5,
              xmax=i+0.5, linestyles='dashed',
              color='black')

plt.xlabel('')
ax.set_xticklabels(ax.get_xticklabels(),rotation = 30)
plt.tight_layout()
fig.show()
fig.savefig('quality_db_homologs.svg')
fig.savefig('quality_db_homologs.png', dpi=600)

