import logging
import pandas as pd
from tqdm import tqdm
from itertools import combinations
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests as mt

logging.warning('# load data')
data = pd.read_table('data/complete_amps_associated_taxonomy.tsv.gz')
gen = pd.read_table('data/general_envo_names.tsv.xz')
metadata = pd.read_table("data/metadata.tsv.xz")

metadata = metadata.merge(on=['microontology',
                              'host_scientific_name',
                              'host_tax_id'],
                          right=gen,
                          how='outer')

metadata = metadata[~metadata.sample_accession.isna()]

metadata = metadata[['sample_accession',
                     'geographic_location',
                     'latitude',
                     'longitude',
                     'general_envo_name',
                     'environment_material',
                     ]]

metadata = metadata.rename({'sample_accession': 'sample'}, axis=1)

logging.warning('# filter and add metadata info')
data = data[data.level == 'species']

logging.warning('# clear duplicates and count nr_AMPs per species per habitat')
amp_rich = data[['source', 'sample', 'amp']]
amp_rich = amp_rich.drop_duplicates()
amp_rich = amp_rich.groupby(['source', 'sample'])
amp_rich = amp_rich.agg('size')
amp_rich = amp_rich.reset_index()
amp_rich = amp_rich.merge(on='sample', right=metadata)
amp_rich = amp_rich.rename({0: 'nr_amp'}, axis=1)

logging.warning('# fix habitat name')
amp_rich.general_envo_name = amp_rich.general_envo_name.replace('human saliva', 'human mouth')

logging.warning('# select species present in at >=10 samples in a habitat')
k = amp_rich.groupby(['source', 'general_envo_name']).agg('size')
k = k[k >= 10]

logging.warning('# select species present in more than 1 habitat')
k = k.reset_index().groupby('source').apply(lambda x: x.general_envo_name.tolist())
k = k[k.apply(lambda x: len(x)) > 1]
k = k.reset_index()

logging.warning('# test difference in the nr_AMPs among samples of the same species through different environments')
test = []
for _, s, env in k.itertuples():
    combs = combinations(env, 2)
    for i, j in combs:
        print(f'{s}\t{i}\t{j}')
        n1 = amp_rich.loc[(amp_rich.source == s) & (amp_rich.general_envo_name == i), 'nr_amp']
        n2 = amp_rich.loc[(amp_rich.source == s) & (amp_rich.general_envo_name == j), 'nr_amp']
        u, p = mannwhitneyu(n1, n2)
        test.append((s, i, len(n1),
                     n1.mean(), n1.std(),
                     j, len(n2), n2.mean(),
                     n2.std(), p))

logging.warning('# format results')
test = pd.DataFrame(test, columns=['species',
                                   'env1',
                                   'n1_samples',
                                   'amp1_avg',
                                   'amp1_std',
                                   'env2',
                                   'n2_samples',
                                   'amp2_avg',
                                   'amp2_std',
                                   'p_value'])

logging.warning('# adjust p-values')
_, test['p_adj'], _, _ = mt(test['p_value'], method='hs')

logging.warning('# export')
test.to_csv('species_amp_richness_crossenvironment.tsv.gz',
            sep='\t',
            header=True,
            index=None)

