import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sp_class = pd.read_table('results_per_outcome.tsv')
t = sp_class[['species', '#_amps', 'outcome']].sort_values(by='outcome')

for val in t.groupby('outcome'):
    val[1]['#_amps'].plot.hist(bins=12, label=val[0], alpha=0.25)

plt.legend()
plt.xlabel('AMPs per genome')
plt.savefig('hist_dist.svg')
plt.close()

for val in t.groupby('outcome'):
    sns.kdeplot(data=val[1], x='#_amps', label=val[0], alpha=0.5)

plt.xlim(0,10)
plt.xlabel('AMPs per genome')
plt.legend()
plt.savefig('density_dist.svg')
plt.close()

for val in t.groupby('species'):
    for j in val[1].groupby('outcome'):
        sns.kdeplot(data=j[1], x='#_amps', label=j[0], alpha=0.5)
    plt.xlim(0,10)
    plt.xlabel('AMPs per genome')
    plt.legend()
    plt.savefig(f'density_{val[0]}_dist.svg')
    plt.close()


t.outcome = t.outcome.apply(lambda x: x.split(' ')[0])
t.outcome = t.outcome.replace('species', 'lost')

