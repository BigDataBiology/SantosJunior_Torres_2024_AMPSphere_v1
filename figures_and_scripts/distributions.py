import pandas as pd
import plotly.express as px
from plotly.graph_objects import Figure
import livingTree as lt
import datatable as dt


def main():
    print('Loading tables...')
    data_prefix = '/home1/huichong/SyncAWS/AMPSphereMetadata/'
    metadata = dt.fread(data_prefix + '/Metadata.tsv', sep='\t', header=False).to_pandas()
    metadata.columns = 'gmsc, amp, sample, microbial_source, specI, is_metagenomic, geographic_location, ' \
                       'latitude, longitude, general_envo_name, environment_material'.split(', ')
    print(metadata)
    print('Start processing...')
    color_map = {}
    metadata['latitude'] = metadata['latitude'].round(0)
    metadata['longitude'] = metadata['longitude'].round(0)
    # metadata['habitat_type'] = pd.Categorical(metadata['microontology'].apply(lambda x: x.split(':')[0]))
    # TODO add legend for the map form number of amps to the bubble size.
    # TODO color consistency across the geo distribution map and habitat distribution.
    #metadata['color'] = metadata['habitat_type'].map(color_map)

    def simplify(data, type='microbial_source'):
        data = data.sort_values(by='size', ascending=False)
        top_9 = data[0:9]
        return pd.DataFrame(top_9.values.tolist(), columns=data.columns).\
            append({type: 'others', 'size': data.loc[9:, 'size'].sum()}, ignore_index=True)

    data = dict(
        geo=metadata[['amp', 'latitude', 'longitude']].
            groupby(['latitude', 'longitude'], as_index=False, observed=True).size(),
        habitats=simplify(metadata[['amp', 'general_envo_name']].groupby('general_envo_name', as_index=False).size(),
                          'general_envo_name'),
        origins=simplify(metadata[['amp', 'microbial_source']].groupby('microbial_source', as_index=False).size(),
                         'microbial_source')
    )

    print('Processing geo data...')
    names = {'latitude': 'lat', 'longitude': 'lon', 'amp': 'size'}
    data['geo'].rename(columns=names, inplace=True)
    # TODO increase the size of bubble
    # data['geo']['size'] = data['geo']['size'] * 5
    print(data['geo'])

    print('Generating plots...')
    figures = dict(
        geo=GeoDistribution(data['geo'], lat='lat', lon='lon', size='size', color='habitat_type'),
        habitats=px.bar(data['habitats'], y='general_envo_name', x='size', orientation='h', color_discrete_sequence=['#1b9e77']),
        origins=px.bar(data['origins'], y='microbial_source', x='size', orientation='h', color_discrete_sequence=['#d95f02'])
    )

    config = {
        'toImageButtonOptions': {
            'format': 'svg',  # one of png, svg, jpeg, webp
            'filename': 'custom_image',
            'height': 500,
            'width': 700,
            'scale': 3  # Multiply title/legend/axis/canvas sizes by this factor
        }
    }

    print('Saving plots...')
    figures['geo'].write_html('geoDistribution.html', config=config)
    figures['habitats'].write_html('habitatDistribution.html', config=config)
    figures['origins'].write_html('originDistribution.html', config=config)
    # figures['geo'].write_image('geoDistribution.svg')
    # figures['habitats'].write_image('habitatDistribution.svg')
    # figures['origins'].write_image('originDistribution.svg')
 

def GeoDistribution(data, lat, lon, size, color) -> Figure:
    fig = px.scatter_geo(data,
                         lat=lat,
                         lon=lon,
                         #color=color,  TODO fixme
                         locationmode='ISO-3',
                         size=size)
    fig.update_layout(
        geo=dict(
            scope='world',
            landcolor='rgb(217, 217, 217)')
    )
    return fig


if __name__ == '__main__':
    main()
