import pandas as pd
from collections import Counter


# define function to classify the lowest taxonomical level associated to a given AMP
def selector(i: tuple):
    if i[1]['level'].isin(['species']).sum() == 1: return f'{i[0]}\tspecies'
    if i[1]['level'].isin(['genus']).sum() == 1: return f'{i[0]}\tgenus'
    if i[1]['level'].isin(['family']).sum() == 1: return f'{i[0]}\tfamily'
    if i[1]['level'].isin(['order']).sum() == 1: return f'{i[0]}\torder'
    if i[1]['level'].isin(['class']).sum() == 1: return f'{i[0]}\tclass'
    if i[1]['level'].isin(['phylum']).sum() == 1: return f'{i[0]}\tphylum'
    if i[1]['level'].isin(['superkingdom']).sum() == 1: return f'{i[0]}\tsuperkingdom'
 

# file available in > ubuntu@aws.big-data-biology.org:/share/work/Celio/files_for_figures/genes_origins
data = pd.read_table('complete_amps_associated_taxonomy.tsv.gz', sep='\t', header='infer')

# filter dataframe and drop duplicates
data[['gmsc', 'amp', 'sample', 'taxid', 'level', 'name', 'specI']]
dd = data[['amp', 'level', 'name', 'specI']].drop_duplicates()

print('# number of AMPs annotated at some level')
print(dd[~dd['level'].isna()]['amp'].drop_duplicates())

# filter AMPs not associated to any origin
ds = dd[~dd['level'].isna()][dd['level'] != 'no rank'][['amp', 'level']]
ds = ds.drop_duplicates()
ds = ds.sort_values(by=['amp', 'level'])
ds = ds.reset_index(drop=True)

# classify AMPs according lowest known levels   
with open('result_tax.txt', 'w') as ofile:
    ofile.write('amp\ttaxonomy\n')
    for i in ds.groupby('amp'):
        ofile.write(f'{selector(i)}\n')
        
# save results as a table
da = pd.read_table('result_tax.txt', sep='\t', header='infer')
pd.DataFrame.from_dict(Counter(da.taxonomy), orient='index')
selected_group = da[da.taxonomy.isin(['species', 'genus'])]['amp']

# convert dataframe of AMPs with annotated origins
sp = data[(data['amp'].isin(selected_group)) & (data['level'].isin(['species', 'genus']))]
sp = sp[['amp', 'taxid', 'level', 'name']]
sp = sp.sort_values(by=['amp', 'taxid', 'name'])
sp = sp.drop_duplicates()
sp['fixed'] = [x.split(' ')[0] for x in sp.name.values]
sp.to_csv('taxonomy_annotation.tsv', sep='\t', header=True, index=None)

## finding amPs annotated to more than 1 genus:
sp_next = sp[['amp', 'fixed']].drop_duplicates()
Counter(Counter(sp_next.amp).values())
# Counter({1: 507470, 2: 53609, 3: 6437, 4: 1594, 5: 539, 6: 219, 7: 99, 8: 80, 9: 35, 10: 24, 11: 17, 12: 17, 13: 11, 15: 8, 14: 6, 16: 5, 23: 4, 17: 4, 26: 2, 22: 2, 21: 2, 24: 1, 18: 1, 20: 1})

## finding 25 most common genera:
ps = pd.DataFrame.from_dict(Counter(sp_next.fixed), orient='index').sort_values(0) * 100 / 863498
ps.tail(25).plot.bar(rot=30, legend=False)
plt.ylabel('% of AMPSphere candidates')
plt.xlabel('Top 25 genera in AMPSphere')
plt.tight_layout()
plt.savefig('top25genera.svg')

## finding 10 most common genera:
ps.tail(10).plot.bar(rot=30, legend=False)
plt.ylabel('% of AMPSphere candidates')
plt.xlabel('Top 10 genera in AMPSphere')
plt.tight_layout()
plt.savefig('top10genera.svg')

