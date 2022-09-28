def plot_mundi_map():
    '''
    It plots the mundi map of all metagenome samples
    :return: A map with all samples containing AMPs
             displayed as blue dots
    '''
    import pandas as pd
    import geopandas as gpd
    import matplotlib.pyplot as plt

    print('Import data')
    input_file = 'data/gmsc_amp_genes_envohr_source.tsv.gz'

    data = pd.read_table(input_file,
                         sep='\t',
                         header='infer')

    print('# filter coordinates and eliminate redundancy')
    df = data[['latitude', 'longitude']]
    df = df.dropna()
    df = df.drop_duplicates()
    df = df.astype('float')
    
    print('# generate background')
    countries = gpd.read_file(
        gpd.datasets.get_path("naturalearth_lowres"))

    print('# initialize an axis')
    fig, ax = plt.subplots(figsize=(8, 6))

    print('# plot map on axis')
    countries.plot(color="lightgrey", ax=ax)

    print('# plot points')
    df.plot(x="longitude", y="latitude", kind="scatter",
            colormap="YlOrRd", s=0.5, ax=ax)

    print('# savefigure')
    plt.savefig('figure_1a_metagenomes_samples_distribution.png',
                dpi=1200)
    plt.close()

