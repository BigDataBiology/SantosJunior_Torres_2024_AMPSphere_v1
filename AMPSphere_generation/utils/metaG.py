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
    
    corrlist = f'{analysis_folder}/AMPsphere_GMSC_correspondence.tsv.gz'
    corrlist = pd.read_table(corrlist, sep='\t', header='infer')

    # generating correspondence list of access and genes
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
    for record in pd.read_table(f'{data_folder}/shortened_renaming.txt.gz',
                                sep='\t',
                                header=None,
                                chunksize=2_500_000):
        record = record.merge(on=0, right=dictresource)
        for f in record.itertuples():
            f = format_geneinfo(f)
            if len(f) == 13:
                df = pd.DataFrame(f).T
                newdf = pd.concat([newdf, df])

    newdf.columns = columns
    newdf.to_csv(f'{analysis_folder}/AMPsphere_metaG_annotation.tsv.gz',
                 sep='\t', header=True, index=None)


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
    species_out = f'{analysis_folder}/AMPSphere_v.2021-03.species.tsv.gz'
    origin_file = f'{analysis_folder}/AMPSphere_v.2021-03.origin_samples.tsv.gz'

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
    metagenomes = progenomes.reset_index()
    metagenomes = progenomes.rename({'index': 'accession'},
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
    
    # inputs
    corrlist = f'{analysis_folder}/AMPsphere_GMSC_correspondence.tsv.gz' 
    metagenomes = f'{analysis_folder}/AMPsphere_metaG_annotation.tsv.gz'
    metadata = f'{data_folder}/metadata.tsv'
    # outputs
    host_out = f'{analysis_folder}/AMPSphere_v.2021-03.hosts.tsv.gz'
    location_out = f'{analysis_folder}/AMPSphere_v.2021-03.locations.tsv.gz'
    microont_out = f'{analysis_folder}/AMPSphere_v.2021-03.microontology.tsv.gz'

    


def metag():
    import os
    
    data_folder='data/'
    analysis_folder='analysis/'
    
    for f in [data_folder, analysis_folder]:
        os.makedirs(f, exist_ok=True)

    print('Filter original GMSC large file')    
#    shorten_input(data_folder, analysis_folder)
    print('Work on the filtered subset')
    anno_metaG(data_folder, analysis_folder)
    print('Generating file with AMP origins')
    origins(data_folder, analysis_folder)
    print('Association of AMPs and metadata')
    assoc_metadata(data_folder, analysis_folder)

