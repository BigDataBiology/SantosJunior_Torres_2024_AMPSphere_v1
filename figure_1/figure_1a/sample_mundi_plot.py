import pandas as pd
import geopandas as gpd

# file in > ubuntu@aws.big-data-biology.org:/share/work/Celio/files_for_figures/mundi_map/
data = pd.read_table('gmsc_amp_genes_envohr_source.tsv',
                     sep='\t',
                     header='infer')

# filter coordinates and eliminate redundancy
df = data[['latitude', 'longitude']]
df = df.dropna()
df = df.drop_duplicates()

# generate background
countries = gpd.read_file(
               gpd.datasets.get_path("naturalearth_lowres"))

# initialize an axis
fig, ax = plt.subplots(figsize=(8,6))

# plot map on axis
countries.plot(color="lightgrey", ax=ax)

# plot points
df.plot(x="longitude", y="latitude", kind="scatter", 
        colormap="YlOrRd", s=0.5, ax=ax)

# savefigure
plt.savefig('metagenomes.svg')
plt.close()
