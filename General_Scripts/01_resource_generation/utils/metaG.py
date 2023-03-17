def shorten_input(data_folder, analysis_folder):
    '''
    Filter large input containing all small genes predicted
    from metaG and leave just those predicted in AMPSphere 
    '''
    import pandas as pd

    ampsphere = f'{analysis_folder}/AMPsphere_GMSC_correspondence.tsv.gz'
    infile = f'{data_folder}/GMSC10.metag_smorfs.rename.txt.xz'
    ofile = f'{data_folder}/shortened_renaming.txt.gz'

    # retrieve genes
    data = pd.read_table(ampsphere, sep='\t', header='infer')
    gmsc_genes = []
    data['genes'] = [x.split(',') for x in data['genes']]
    for x in data['genes']: gmsc_genes += x
    gmsc_genes = set(gmsc_genes)
    
    for record in pd.read_table(infile,
                                sep='\t',
                                header='infer', 
                                chunksize=10_000_000):
        record = record[record['#GMSC_id'].isin(gmsc_genes)]
        record.to_csv(ofile,
                      mode='at',
                      sep='\t',
                      header=False,
                      index=None)


def format_geneinfo(prodigal_out):
    '''
    Takes the line from prodigal gene prediction and 
    returns a set of variables

    :Input:
    
    prodigal_out: tuple, 3 fields
                  from pandas itertuples
                  
    :Output:
    
    list, 13 fields 
    '''    
    gmsc, sample = prodigal_out[1], prodigal_out[2]
    contig = '_'.join(prodigal_out[3].split('_')[0:2])
    start, stop = prodigal_out[5], prodigal_out[7]
    strand = prodigal_out[9]
    amp = prodigal_out[12]
    prodigal_out = prodigal_out[11].split(';')
    fid = prodigal_out[0].replace('ID=', '')
    partial = prodigal_out[1].replace('partial=', '')
    start_type = prodigal_out[2].replace('start_type=', '')
    rbs_motif = prodigal_out[3].replace('rbs_motif=', '')
    rbs_spacer = prodigal_out[4].replace('rbs_spacer=', '')
    gc_cont = prodigal_out[5].replace('gc_cont=', '')

    return [amp, gmsc, sample,
            contig, start, stop, strand, 
            fid, partial, start_type,
            rbs_motif, rbs_spacer,
            float(gc_cont)]
        

def anno_metaG(data_folder, analysis_folder):
    '''
    Fix the format for the shortened metaG info input
    '''
    import pandas as pd

    # generating correspondence list of access and genes    
    corrlist = f'{analysis_folder}/AMPsphere_GMSC_correspondence.tsv.gz'
    corrlist = pd.read_table(corrlist, sep='\t', header='infer')
    corrlist = corrlist[['accession', 'genes']]
    corrlist = corrlist.set_index('accession')
    corrlist['genes'] = [x.split(',') for x in corrlist.genes]

    dictresource = dict()
    for f in corrlist.itertuples():
        for x in f[1]: dictresource[x] = f[0]

    dictresource = pd.DataFrame.from_dict(dictresource, orient='index')
    dictresource = dictresource.reset_index()
    dictresource = dictresource.rename({'index': 0, 0: 1}, axis=1)
    
    columns = ['accession', 'gmsc', 'sample', 'contig', 'start', 'stop', 'strand', 
               'fid', 'partial', 'start_type', 'rbs_motif',
               'rbs_spacer', 'gc_cont']

    newdf = pd.DataFrame()
    infile = f'{data_folder}/shortened_renaming.txt.gz'
    for idx, record in enumerate(pd.read_table(infile,
                                               sep='\t',
                                               header=None,
                                               chunksize=500_000)):
        record = record.merge(on=0, right=dictresource)
        record = [format_geneinfo(f) for f in record.itertuples()]
        df = pd.DataFrame(record, columns=columns)
        newdf = pd.concat([newdf, df])

    ofile = f'{analysis_folder}/AMPsphere_metaG_annotation.tsv.gz'
    
    newdf = newdf.sort_values(by=['accession',
                                  'gmsc',
                                  'sample'],
                              ascending=[True,
                                         True,
                                         True])
                                         
    newdf.to_csv(ofile,
                 sep='\t',
                 header=True,
                 index=None)


def origins(data_folder, analysis_folder):
    '''
    Annotate AMPs origins discriminating if it comes from
    metagenomes, progenomes, and their genes
    '''
    import pandas as pd
    
    # inputs
    corrlist = f'{analysis_folder}/AMPsphere_GMSC_correspondence.tsv.gz' 
    progenomes = f'{analysis_folder}/AMPsphere_proGenomes_correspondence.tsv.gz'
    metagenomes = f'{analysis_folder}/AMPsphere_metaG_annotation.tsv.gz'
    # outputs
    species_out = f'{analysis_folder}/AMPSphere_v.2022-03.species.tsv.gz'
    origin_file = f'{analysis_folder}/AMPSphere_v.2022-03.origin_samples.tsv.gz'

    # processing gmsc genes
    gmsc = pd.read_table(corrlist, sep='\t', header='infer')

    # processing progenomes
    data = pd.read_table(progenomes, sep='\t', header='infer')
    data[['accession',
          'genome',
          'specI']].to_csv(species_out, 
                           sep='\t', header=True, index=None)

    progenomes = dict()
    for record in data.groupby('accession'):
        progenomes[record[0]] = ','.join(record[1].genome.tolist())

    progenomes = pd.DataFrame.from_dict(progenomes,
                                        orient='index',
                                        columns=['progenomes'])
                                        
    progenomes = progenomes.reset_index()
    progenomes = progenomes.rename({'index': 'accession'},
                                   axis=1)

    # processing metagenomes
    data = pd.read_table(metagenomes, sep='\t', header='infer')
    metagenomes = dict()
    for record in data.groupby('accession'):
        metagenomes[record[0]] = ','.join(record[1]['sample'].tolist())
    
    metagenomes = pd.DataFrame.from_dict(metagenomes,
                                         orient='index',
                                         columns=['metagenomes'])
                                        
    metagenomes = metagenomes.reset_index()
    metagenomes = metagenomes.rename({'index': 'accession'},
                                     axis=1)
    
    # merging
    gmsc = gmsc.merge(on='accession', right=progenomes, how='outer')
    gmsc = gmsc.merge(on='accession', right=metagenomes, how='outer')
    
    gmsc = gmsc[['accession', 'genes', 'progenomes', 'metagenomes']]
    
    gmsc.to_csv(origin_file, 
                sep='\t',
                header=True,
                index=None)


def assoc_metadata(data_folder, analysis_folder):
    '''
    Associate metadata to the AMPSphere info
    '''
    import pandas as pd
    from collections import Counter
    
    # inputs
    metagenomes = f'{analysis_folder}/AMPsphere_metaG_annotation.tsv.gz'
    metadata = f'{data_folder}/metadata.tsv.xz'

    # outputs
    host_out = f'{analysis_folder}/AMPSphere_v.2022-03.hosts.tsv.gz'
    location_out = f'{analysis_folder}/AMPSphere_v.2022-03.locations.tsv.gz'
    microont_out = f'{analysis_folder}/AMPSphere_v.2022-03.microontology.tsv.gz'

    # load metagenomes
    metagenomes = pd.read_table(metagenomes, sep='\t', header='infer')
    metagenomes = metagenomes[['accession', 'sample']]
    
    # load metadata
    metadata = pd.read_table(metadata, sep='\t', header='infer')
    metadata = metadata.rename({'sample_accession': 'sample'}, axis=1)
    metadata = metadata[['sample',
                         'microontology',
                         'geographic_location',
                         'latitude',
                         'longitude',
                         'environment_material',
                         'host_common_name',
                         'host_scientific_name',
                         'host_tax_id']]

    # merge
    metagenomes = metagenomes.merge(on='sample', right=metadata)
    
    # microontology
    micro = metagenomes[['accession',
                         'microontology']]
    
    micro = micro.pivot_table(index=['accession',
                                     'microontology'],
                              aggfunc='size')

    micro = micro.reset_index()
    micro = micro.rename({0: 'counts'},
                         axis=1)                              
                         
    micro.to_csv(microont_out,
                 sep='\t',
                 header=True,
                 index=None)    

    # geoloc
    geoloc = metagenomes[['accession',
                          'geographic_location']]
    
    geoloc = geoloc.pivot_table(index=['accession',
                                       'geographic_location'],
                                aggfunc='size')
                                
    geoloc = geoloc.reset_index()
    geoloc = geoloc.rename({0: 'counts'},
                         axis=1)                              

    geoloc.to_csv(location_out,
                  sep='\t',
                  header=True,
                  index=None)    

    # host_out
    host = metagenomes[['accession',
                        'host_common_name',
                        'host_scientific_name',
                        'host_tax_id']]

    host = host.dropna()                    
    host = host.pivot_table(index=['accession',
                                   'host_common_name',
                                   'host_scientific_name',
                                   'host_tax_id'],
                            aggfunc='size')
                            
    host = host.reset_index()
    host = host.rename({0: 'counts'},
                       axis=1)                              

    host.to_csv(host_out,
                sep='\t',
                header=True,
                index=None)


def metag():
    import os
    
    data_folder='data/'
    analysis_folder='analysis/'
    
    for f in [data_folder, analysis_folder]:
        os.makedirs(f, exist_ok=True)

    print('Filter original GMSC large file')    
    shorten_input(data_folder, analysis_folder)
    print('Work on the filtered subset')
    anno_metaG(data_folder, analysis_folder)
    print('Generating file with AMP origins')
    origins(data_folder, analysis_folder)
    print('Association of AMPs and metadata')
    assoc_metadata(data_folder, analysis_folder)

