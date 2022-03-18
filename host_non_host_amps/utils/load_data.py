def load_data():
    ## loading data
    import pandas as pd
    print('... import taxonomy')
    tax = pd.read_table('data/taxonomy_annotation.tsv')
    tax = tax[['amp', 'fixed']]
    tax = tax.drop_duplicates()
    print('... iterate over groups')
    tax = tax.groupby('amp').agg(lambda x: set(x))
    tax = tax.reset_index()
    print('... load families')
    spheres = pd.read_table('data/SPHERE_v.2022-03.levels_assessment.tsv.gz')
    spheres = spheres[['AMP accession', 'SPHERE_fam level III']]
    spheres = spheres.rename({'AMP accession': 'amp',
                              'SPHERE_fam level III': 'family'},
                             axis=1)
    return tax, spheres


def load_general_envo(tax):
    ## stating environmental breaking up:
    from .infolists import mammal, nonmammal, plant, environmental
    import pandas as pd
    print('... import GMSC metadata info')
    data = pd.read_table('data/gmsc_amp_genes_envohr_source.tsv.gz')
    print('... eliminate AMP rows belonging to ProGenomes')
    data = data[(~data['general_envo_name'].isin(['N.A.','*']))]
    data = data[(~data['general_envo_name'].isna())]  
    print('... reduces size of table')
    data = data[['amp', 'sample', 'general_envo_name']]  
    data = data.drop_duplicates()  
    data = data.merge(on='amp', right=tax)
    print('... environmental filter')
    data['general_envo_name'] = data['general_envo_name'].replace(mammal, 'mammal')
    data['general_envo_name'] = data['general_envo_name'].replace(nonmammal, 'non-mammal')
    data['general_envo_name'] = data['general_envo_name'].replace(plant, 'plant')
    data['general_envo_name'] = data['general_envo_name'].replace(environmental, 'environmental')
    print('... fixing column')
    data.fixed = [str(x) for x in data.fixed]
    return data


def input_info():
    tax, spheres = load_data()
    data = load_general_envo(tax)
    return data, spheres
    
