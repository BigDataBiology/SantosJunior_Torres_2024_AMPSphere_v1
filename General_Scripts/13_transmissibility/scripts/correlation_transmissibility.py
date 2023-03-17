import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import spearmanr


# getting the columns to correlate in the final test
cols = ['SGB_mother-offspring_transmissibility_until1y',
        'SGB_mother-offspring_transmissibility_until3y',
        'SGB_mother-offspring_transmissibility_until18y',
        'SGB_household_transmissibility',
        'SGB_intradataset_transmissibility', 
        'SGB_mother-infant_transmissibility']


def convdict(x: str):
    '''
    convdict is a function to convert a taxonomy lineage
    string into a dictionary returning the corresponding
    species or a None
    
    input:
    
    taxonomy lineage string with fields separated by | 
    with the different taxonomy levels separated by __
    
    output:
    
    string corresponding to the species or None in case
    it is not available
    '''
    f = dict()
    x = x.split('|')
    for j in x:
        j = j.split('__')
        f[j[0]] = j[1]
    return f.get('s', None)


def corrs(df, output:str):
    '''
    corrs - Calculates and plot the correlations between 
    species density of AMPs and their transmissibility
    
    input:
        df = dataframe of pandas containing different columns
        with transmissibility measures and one called density
        containing number of AMP genes per assembled base pairs
        also needs to contain a column named species, which is
        used as reference in correlation
        
        output = string used to save the data
    
    output:
        corrs = dataframe object with correlation measures
    '''
    corrs = []   
    for j in cols:
        if j in df.columns:
            r, p = spearmanr(df['density'],
                             df[j],
                             nan_policy='omit')
            L = len(df[j].dropna())
            corrs.append((j, r, p, L))
    corrs = pd.DataFrame(corrs,
                         columns=['Transmissibility', 'SpearmanR',
                                  'p-value', 'N'])
    for j in cols:
        if j in df.columns:
            x = df[['density', j]].dropna()
            x = x['density']
            y = df[j].dropna()
            sns.scatterplot(x=x, y=y)
            plt.xscale('log')
            plt.xlabel("AMP density per Mbp")
            plt.yscale('log')
            plt.ylabel("Transmissibility")
            plt.tight_layout()
            plt.savefig(f'{output}_{j}.svg')
            plt.close()
    df.to_csv(f'results_table_{output}.tsv',
              sep='\t',
              header=True,
              index=True)           
    corrs.to_csv(f'results_corr_{output}.tsv',
                 sep='\t',
                 header=True,
                 index=None)
    return corrs         


# getting species AMP density
dens = pd.read_table('species_amp_density_per_sample.tsv.gz')
dens = dens[dens.general_envo_name.str.contains('human')]  # filter only human-associated envo
dens = dens.groupby('name').apply(lambda x: x[['amp_genes', 'nbps']].sum(axis=0))  # sum values for the same species
dens['density'] = 1_000_000 * dens.amp_genes / dens.nbps  # generates density per Mbp
dens = dens.reset_index(drop=False)
dens.rename({'name': 'species'}, axis=1, inplace=True)  # fix column name to merging it later

# getting species from the paper
#s4 - description of taxonomy for stool species
s4 = pd.read_excel('41586_2022_5620_MOESM3_ESM.xlsx',
                   sheet_name='Table S4')
s4['species'] = s4.Taxonomy.apply(lambda x: convdict(x))
s4 = s4[['SGB', 'species']]
s4 = s4[s4.species != None]
s4.species = s4.species.apply(lambda x: ' '.join(x.split('_')[0:2]))

# getting transmissibility measures for gut samples
s9 = pd.read_excel('41586_2022_5620_MOESM3_ESM.xlsx', sheet_name='Table S9')
s4 = s4.merge(on='SGB', right=s9)
s4 = s4.merge(on='species', right=dens)
corrs(s4, 'gut')

## Results for gut:
#                    Transmissibility  SpearmanR   p-value    N
#      SGB_household_transmissibility  -0.129400  0.201767   99
#   SGB_intradataset_transmissibility  -0.049839  0.596837  115
#  SGB_mother-infant_transmissibility  -0.440864  0.002436   45

#s5 - description of taxonomy for saliva species
s5 = pd.read_excel('41586_2022_5620_MOESM3_ESM.xlsx',
                   sheet_name='Table S5')
s5['species'] = s5.Taxonomy.apply(lambda x: convdict(x))
s5 = s5[['SGB', 'species']]
s5 = s5[s5.species != None]
s5.species = s5.species.apply(lambda x: ' '.join(x.split('_')[0:2]))

# getting transmissibility measures for oral samples
s30 = pd.read_excel('41586_2022_5620_MOESM3_ESM.xlsx', sheet_name='Table S30')
s5 = s5.merge(on='SGB', right=s30)
s5 = s5.merge(on='species', right=dens)
corrs(s5, 'oral')

## Results for oral:
#                                Transmissibility  SpearmanR   p-value   N
#   SGB_mother-offspring_transmissibility_until1y  -0.442549  0.002072  46
#   SGB_mother-offspring_transmissibility_until3y  -0.308567  0.005662  79
#  SGB_mother-offspring_transmissibility_until18y  -0.260272  0.010856  95
#                  SGB_household_transmissibility  -0.211325  0.039806  95
#               SGB_intradataset_transmissibility   0.159904  0.121647  95

