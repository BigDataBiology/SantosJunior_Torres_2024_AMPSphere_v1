def create_GMGCgenes():
    '''
    Creates a list of GMSC genes
    '''
    import pandas as pd
    
    data = pd.read_table('data/AMPSphere_v.2021-03.origin_samples.tsv.gz', sep='\t', header='infer')

    gmscgenes = []
    for i in data['GMSC accession']:
        gmscgenes += i.split(',')    

    gmscgenes.sort()
    with open('data/GMGCgenes', 'w') as ofile:
       for i in gmscgenes:
           ofile.write(f'{i}\n')

    return gmscgenes


def filtergenes(gmscgenes):
    '''
    Filter gene sequences by code
    and retrieve their coding info
    '''
    import lzma
    import pandas as pd

    fdf = pd.DataFrame()
    for j, dv in enumerate(pd.read_table('data/GMSC10.metag_smorfs.rename.txt.xz',
                                         sep='\t',
                                         header='infer',
                                         chunksize=10_000_000)):
        df = dv[dv['#GMSC_id'].isin(gmscgenes)]
        fdf = pd.concat([fdf, df], ignore_index=True)
        print(f'... processed chunk {j}')

    fdf.to_csv('data/gmsc_genes_smorfs.txt', sep='\t', header=True, index=None)

    return fdf


def header_work(header: str):
    '''
    Break the header using format patterns into different 
    classes of variables
    '''
    header = header.replace(' # ', ';')
    header = header.split(';')
    c[0] = '_'.join(c[0].split('_')[:-1])
    c[4] = c[4].replace('ID=', '')
    c[5] = c[5].replace('partial=', '')
    c[6] = c[6].replace('start_type=', ''))
    c[7] = c[7].replace('rbs_motif=', '')
    c[8] = c[8].replace('rbs_spacer=', '')
    c[9] = c[9].replace('gc_cont=', '')
    
    # return contig, start, stop, strand, ID, partial, 
    #        start_type, rbs_motif, rbs_spacer, gc_cont
    return (c[0], c[1], c[2], c[3],
            c[4], c[5], c[6], c[7],
            c[8], c[9])


def genes2amp():
    '''
    Creates a dictionary relating 
    the gmsc gene and the AMP in 
    AMPSphere
    '''
    from Bio import SeqIO
    import lzma

    genes_to_amp = dict()
    ampsphere = lzma.open('data/AMPSphere_v.2021-03.fna.xz', 'rt', encoding='utf-8')
    for record in SeqIO.parse(ampsphere, 'fasta'):
        genes_to_amp[record.description.split()[0]] = record.description.split()[1]

    ampsphere.close()
    return genes_to_amp


def format_smorfs(fdf):
    '''
    Take the gmsc genes info and format the dataframe
    into an specific set of columns
    '''
    import pandas as pd
    import numpy as np
    from tqdm import tqdm

    print('... preparing arrays')
    contig, start, stop, strand = [],[],[],[]
    ID, partial, start_type, rbs_motif = [],[],[],[]
    rbs_spacer, gc_cont = [],[]

    print('... separating and classifying fdf')
    with tqdm(total=len(fdf),
              desc="Processed pairs",
              bar_format="{l_bar}{bar} [ time left: {remaining} ]") as pbar:
        for c in fdf.original_header.tolist():
            c = header_work(c)
            contig.append(c[0])
            start.append(c[1])
            stop.append(c[2])
            strand.append(c[3])
            ID.append(c[4])
            partial.append(c[5])
            start_type.append(c[6])
            rbs_motif.append(c[7])
            rbs_spacer.append(c[8])
            gc_cont.append(c[9])
            pbar.update(1)
        
    print('... creating aditional classes')
    gmsc = fdf['#GMSC_id'].tolist()
    samples = fdf['sample'].tolist()
    
    print('... associating amp names')
    genes_to_amp = genes2amp()
    
    amps = [genes_to_amp.get(gene,
                             'NA') for gene in gmsc]

    print('... forming fdf frame')

    del fdf, genes_to_amp

    df = pd.dataFrame([gmsc, amps, samples, contig,
                       start, stop, strand,
                       ID, partial, start_type,
                       rbs_motif, rbs_spacer, gc_cont],
                      columns=['gmsc', 'amp', 'sample', 'contig',
                               'start', 'stop', 'strand',
                               'ID', 'partial', 'start_type',
                               'rbs_motif', 'rbs_spacer', 'gc_cont'])
    
    print('... fdf exporting')
    df.to_csv('data/gmsc_metag_amps.tsv.gz',
              sep='\t',
              header=True,
              index=None)

    return df
    

def filtermmseqs2_amp_contig(df):
    '''
    Filter the taxonomy annotated with mmseqs2 and GTDB
    The results of the mmseqs2 run are processed using
    the contig names to generate a smaller file
    that will be used later to annotate the AMPs
    '''
    import pandas as pd

    print('... filtering na')
    df = df[df['amp'] != 'NA']

    print('... cleaning tuples')
    df = df[['sample', 'contig']]

    df = df.apply(tuple, axis=1)

    print('... proceding to the filtering')
    fdf = pd.DataFrame()

    for j, dv in enumerate(pd.read_table('data/mmseqs2.lca_taxonomy.full.tsv.xz',
                                         sep='\t',
                                         header='infer',
                                         chunksize=10_000_000)):
        dfd = dv[dv[['sample', 'contig']].apply(tuple, axis=1).isin(df)]
        fdf = pd.concat([fdf, dfd], ignore_index=True)
        print(f'... processed chunk {j}')

    print('... data exporting')
    fdf.to_csv('data/amp_contigs_filtered_mmseqs2.lca_taxonomy.tsv.xz',
               sep='\t',
               header=True,
               index=None)
    
    return fdf 


def add_taxa(df1):
    '''
    Add the taxonomic affiliation to the gene info
    creating a complete table with all required info
    to complete trace back AMP genes origin and structure

    Inputs:

    df1 - filtered results table of taxonomy (dataframe: amps_contigs_...lca_taxonomy.tsv.xz)
    '''
    import pandas as pd

    print('...loading genes')
    df2 = pd.read_table('data/gmsc_metag_amps.tsv.gz', sep='\t', header='infer')
    df2 = df2.dropna()

    print('...merging taxonomies')
    df3 = df2.merge(on=['sample', 'contig'], right=df1)
    df4 = df2[~df2['gmsc'].isin(df3['gmsc'])]

    print('...adding zero taxonomy')
    df5 = pd.concat([df3, df4])
    df5 = df5.sort_values(by='gmsc')

    print('...exporting data')
    df5.to_csv('data/gmsc_meta_taxo.tsv.gz',
               sep='\t',
               header=True,
               index=None)
    

def preprocess():
    print('Generating genes list')
    gmscgenes = create_GMGCgenes()
    print('Filter gene features')
    fdf = filtergenes(gmscgenes)
    print('Format smORFs')
    df = format_smorfs(fdf)
    print('Filter results from mmseqs2')
    taxfilt = filtermmseqs2_amp_contig(df)
    print('Associating taxa to genes')
    add_taxa(taxfilt)
    
