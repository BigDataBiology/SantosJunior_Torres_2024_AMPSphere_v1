import pandas as pd
import plotly.express as px
from plotly.graph_objects import Figure
import livingTree as lt
import datatable as dt


def main():
    print('Loading tables...')
    data_prefix = '/home1/huichong/SyncAWS/AMPSphereMetadata/outputs'
    metadata = dt.fread(data_prefix + '/AMPSphere_metadata.tsv.gz', sep='\t', header=True).to_pandas()
    print(metadata)
    """
                AMPSphere_code                          GMGC          sample  ... host_scientific_name   latitude   longitude
    0        AMP10.000_000  GMSC10.SMORF.000_036_844_556  SAMEA104142073  ...         Homo sapiens  23.127100  113.282800
    1        AMP10.000_000  GMSC10.SMORF.000_036_899_109  SAMEA104142074  ...         Homo sapiens  23.127100  113.282800
    2        AMP10.000_000  GMSC10.SMORF.000_036_944_928  SAMEA104142075  ...         Homo sapiens  23.127100  113.282800
    """
    origins = dt.fread(data_prefix + '/association_table_origins.tsv.gz', sep='\t', header=True).to_pandas()
    print(origins)
    """
             AMPsphere code                          GMSC          sample  taxid                    name
    0        AMP10.000_000  GMSC10.SMORF.000_036_844_556  SAMEA104142073    237             Phocaeicola
    1        AMP10.000_000  GMSC10.SMORF.001_828_692_528    SAMN03955567    237             Phocaeicola
    2        AMP10.000_000  GMSC10.SMORF.001_828_750_506    SAMN03955549    237             Phocaeicola
    """

    print('Start processing...')
    color_map = {}
    metadata['latitude'] = metadata['latitude'].round(0)
    metadata['longitude'] = metadata['longitude'].round(0)
    metadata['habitat_type'] = pd.Categorical(metadata['microontology'].apply(lambda x: x.split(':')[0]))
    #metadata['color'] = metadata['habitat_type'].map(color_map)
    data = dict(
        geo=metadata[['AMPSphere_code', 'latitude', 'longitude', 'habitat_type']].
            groupby(['latitude', 'longitude', 'habitat_type'], as_index=False, observed=True).size(),
        hosts=metadata[['AMPSphere_code', 'host_tax_id', 'host_scientific_name']].
            groupby('host_tax_id', as_index=False).size(),
        habitats=metadata[['AMPSphere_code', 'microontology', 'habitat_type']].
            groupby(['microontology', 'habitat_type'], as_index=False).size(),
        origins=origins[['AMPsphere code', 'taxid', 'name']].
            groupby('taxid', as_index=False).size()
    )

    print('Processing geo data...')
    names = {'latitude': 'lat', 'longitude': 'lon', 'AMPSphere_code': 'size'}
    data['geo'].rename(columns=names, inplace=True)
    print(data['geo'])
    data['habitats'][['l1', 'l2', 'l3', 'l4', 'l5', 'l6']] = data['habitats']['microontology'].str.split(':', expand=True).fillna('Unknown')
    print(data['habitats'])

    print('Assigning lineages to taxa...')
    # Fix id inconsistency.
    data['hosts']['host_tax_id'] = data['hosts']['host_tax_id'].apply(lambda x: x if x != 2116673.0 else 85678.0)
    data['hosts'][['sk', 'k', 'p', 'c', 'o', 'f', 'g', 's']] = \
        pd.DataFrame(lt.LineageTracker(ids=data['hosts']['host_tax_id'].astype(int)).paths_sp,
                     columns=['sk', 'k', 'p', 'c', 'o', 'f', 'g', 's'])
    data['hosts'].fillna('Unknown', inplace=True)

    # data['origins'][['sk', 'k', 'p', 'c', 'o', 'f', 'g', 's']] = \
    #     pd.DataFrame(lt.LineageTracker(ids=data['origins']['taxid'].astype(int)).paths_sp,
    #                  columns=['sk', 'k', 'p', 'c', 'o', 'f', 'g', 's'])
    # data['origins'].fillna('Unknown', inplace=True)
    print('Generating plots...')
    figures = dict(
        geo=GeoDistribution(data['geo'], lat='lat', lon='lon', size='size', color='habitat_type'),
        habitats=sunburstDistribution(data['habitats'], path=['l1', 'l2', 'l3', 'l4', 'l5', 'l6'], value='size'),
        hosts=sunburstDistribution(data['hosts'], path=['sk', 'k', 'p', 'c', 'o', 'f', 'g', 's'], value='size'),
        # origins=sunburstDistribution(data['origins'], path=['sk', 'k', 'p', 'c', 'o', 'f', 'g', 's'], value='size')
    )

    config = {
        'toImageButtonOptions': {
            'format': 'svg',  # one of png, svg, jpeg, webp
            'filename': 'custom_image',
            'height': 500,
            'width': 700,
            'scale': 1  # Multiply title/legend/axis/canvas sizes by this factor
        }
    }

    print('Saving plots...')
    figures['geo'].write_html('geoDistribution.html', config=config)
    figures['habitats'].write_html('habitatDistribution.html', config=config)
    figures['hosts'].write_html('hostDistribution.html', config=config)
    figures['origins'].write_html('originDistribution.html', config=config)


def GeoDistribution(data, lat, lon, size, color) -> Figure:
    fig = px.scatter_geo(data,
                         lat=lat,
                         lon=lon,
                         color=color,
                         locationmode='ISO-3',
                         size=size)
    fig.update_layout(
        geo=dict(
            scope='world',
            landcolor='rgb(217, 217, 217)')
    )
    return fig


def sunburstDistribution(data, path, value) -> Figure:
    fig = px.sunburst(data,
                      path=path,
                      values=value)
    return fig
    pass


if __name__ == '__main__':
    main()
