import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from numpy import mean, std
from scipy.stats import norm
from scipy.stats import kstest
from scipy.stats import shapiro
from scipy.stats import normaltest

df = pd.read_table('data_from_clusters_wo_clones.txt')
df = df[['*' not in x for x in df.NR_proteins_wo_clones]]

core_fams = []
resfams = []
for c in df.Cluster:
    intdf = pd.read_table(f'analysis/{c}_families.tsv.gz')
    x = intdf.loc[intdf.cutoff_prevalence == 95, 'percent of families'].tolist()[0]
    co = intdf.loc[intdf.cutoff_prevalence == 95, 'number of families'].tolist()[0]
    s = intdf.loc[intdf.cutoff_prevalence == 50, 'number of families'].tolist()[0]
    a = intdf.loc[intdf.cutoff_prevalence == 0, 'number of families'].tolist()[0]
    ap = (a-s)*100/a
    sp = (s-co)*100/a
    cop = co*100/a
    core_fams.append(x)
    resfams.append((c, cop, sp, ap))

df = pd.DataFrame(resfams, columns=['cluster', 'core_fams', 'shell_fams', 'accessory_fams'])
df.to_csv('summary_output_core_prots.tsv', sep='\t', header=True, index=None)
print('Generated_table')

df2 = df.set_index('cluster').melt()
sns.boxplot(data=df2, x='variable', y='value', showfliers=False, color='white')
sns.swarmplot(data=df2, x='variable', y='value', s=2.5, alpha=0.5)
plt.xlabel('Conservation across species from ProGenomes2')
plt.xticks([0,1,2], ['Core', 'Shell', 'Accessory'])
plt.ylabel('% of full-length protein families')
plt.savefig('conservation_families.svg')

normality=[]    

_, p = normaltest(core_fams)
normality.append(p)

_, p = shapiro(core_fams)
normality.append(p)

_, p = kstest(core_fams, 'norm')
normality.append(p)

m, s = mean(core_fams), std(core_fams)
zscore = (5.89-m)/s

if any(x > 0.05 for x in normality):
    p = norm().sf(abs(zscore))*2
    print(f'Core families in the proGenomes families follow normal distribution')
    print(f'The z-score for AMP families conservation as core was {zscore:.2f}')
    print(f'That means there is a probability of {p:.2E} of being a random result')
else:
    print(f'Core families in the proGenomes families do not follow normal distribution')
    print(f'The z-score for AMP families conservation as core was {zscore:.2f}')
    pest = 1/(zscore**2)
    print(f'There is an estimate probability by Chebyshev of p < {pest:.2f}')
    p = sum(x<=5.89 for x in core_fams) / len(core_fams)
    print(f'Calculated probability using permutation test was {p:.2f}')

sns.kdeplot(core_fams, clip=(0,100))
plt.axvline(5.89, color='red', linestyle='--')
plt.savefig('core_fams_dist.svg')
plt.close()

sns.boxplot(core_fams, color='white', showfliers=False)
sns.swarmplot(core_fams, color='blue', alpha=0.5, s=3)
sns.stripplot([5.89], color='red', s=5)
plt.savefig('core_fams_dist_box.svg')
plt.close()
