import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import cm
from scipy import stats
from statsmodels.stats import multitest as mt

peptides = pd.read_excel('data/Peptide_nickname_AMPSphere.xlsx', index_col=0)
peptides = peptides[peptides.index.str.startswith('AMP10')].rename(columns={'Active?':'Active'})
peptides['Active'] = peptides['Active'] == 'Yes'
copred = pd.read_table('new_data/AMP_coprediction_AMPSphere.tsv.xz', index_col=0)
quality = pd.read_table('../data_folder/quality_assessment.tsv.xz', index_col=0)

peptides = peptides.join(copred)
peptides = peptides.join(quality)
peptides.eval('hq = Coordinates == "Passed" & Antifam == "Passed" & RNAcode == "Passed"', inplace=True)
peptides.eval('has_meta = metaproteomes == "Passed" | metatranscriptomes == "Passed"', inplace=True)
preds = ['APIN', 'AMPScanner2', 'AMPlify', 'ampir', 'amPEPpy', 'AI4AMP', 'Macrel']

ps = []
_,p = stats.fisher_exact(peptides.groupby(['Active','APIN']).size().values.reshape((2,2)))
ps.append(p)
for p in preds[1:]:
    m_u = (stats.mannwhitneyu(peptides.query('Active')[p], peptides.query('~Active')[p]))
    ps.append(m_u.pvalue)

_, pvals_adjust, _ , _ = mt.multipletests(ps, method='holm-sidak')


fig,ax = plt.subplots()
# convert to long format
# [ Active(bool) | Predictor(str) | value(float) ]
peptides_long = pd.melt(peptides[['Active']+preds], id_vars=['Active'], value_vars=preds)
peptides_long.columns = ['Active', 'Predictor', 'Probability']
peptides_long.Active.replace({True:'Active', False:'Inactive'}, inplace=True)
sns.boxplot(data=peptides_long, x='Predictor', y='Probability', hue='Active', ax=ax, boxprops={'facecolor':'None'}, showfliers=False, legend=False)
s = sns.stripplot(data=peptides_long, x='Predictor', y='Probability', hue='Active', dodge=True, ax=ax, alpha=0.6, palette=cm.Dark2.colors, legend='full')
ax.hlines(0.5, -0.5, len(preds)-.5, linestyles='--', linewidth=1, color='k')
ax.get_legend().set_title(None)
for ix,p in enumerate(pvals_adjust):
    ax.text(ix, 1.05, f'p={p:.02f}', ha='center', va='bottom')
sns.despine(fig, trim=True)
fig.savefig('figures/copredictions.svg')
