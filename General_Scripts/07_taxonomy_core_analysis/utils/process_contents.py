def mapsamples(fin):
    '''
    It maps the sample to its row index. It can be later used to
    select rows in a large xz table. The sample column must be the
    first and the file should be tab separated.
    
    :input:
    - fin     input file opened as text using python-xz package
    
    :output:
    - returns a dictionary mapping the rows and the samples
    '''
    index_sample = dict()
    for idx, row in enumerate(fin):
        if idx > 0:    # assume there is a header
            sample = row.split('\t')[0]
            if sample in index_sample:
                index_sample[sample][1] = idx
            else:
                index_sample[sample] = [idx, idx]
    return index_sample
    

def mapidx(fin, sample_map, ofile):
    '''
    Map the position of each row in the file in bytes
    it is used to random access using python-xz
    
    :input:
    - fin        input file opened as text using python-xz package
    - sample_map input table with sample, start_row, end_row

    :output:
    - returns a dictionary mapping the rows in terms of bytes
    positioning in the file
    '''
    import lzma
    import pandas as pd
    data = pd.read_table(sample_map)
    data['start'] = data['start'] - 1
    data['end'] = data['end'] - 1
    positions = set(data['start']).union(set(data['end']))
    posix = dict()
    idx, row = (0, 0)
    fin.seek(0)
    for l in fin:
        idx += len(l)
        if row in positions: posix[row] = idx
        row += 1
    data['start'] = [posix[x] for x in data['start']]
    data['end'] = [posix[x] for x in data['end']]
    data.to_csv(ofile, sep='\t', header=True, index=None)


def getdf(start: int, end: int, cols: list, fin):
    '''
    Retrieve a dataframe from rows indexes and bytes
    placement in the large file.
    
    :input:
    - start, end      indexes of the first and last rows
                      of interest as bytes in the xz file
    - cols            list of column names in the final
                      dataframe
    - fin             input file opened as text using
                      python-xz package
 
    :output:
    - returns a dataframe in the original format as that
    inputted in the xz file
    '''
    import pandas as pd
    if start == end:
        fin.seek(start)
        f = fin.readline()
        f = f.strip().split('\t')
        f = pd.DataFrame(f, index=cols).T
    else:
        fin.seek(start)
        f = fin.readlines(end-start)
        f = [i.strip().split('\t') for i in f]
        f = pd.DataFrame(f, columns=cols)
    return f


def process_chunk(record, taxrecord):
    '''
    Receives two dataframes, the first containing the number of base pairs of 
    each type and another with the taxonomy of the contigs belonging to that
    sample.
    
    :input:
    - record           dataframe containing the following columns:
                       sample, contig, A, C, T, G
                       
    - taxrecord        dataframe containing the following columns:
                       sample, contig, taxid, level, name, retained,
                       assigned, agreement, support
    
    :output:
    - dataframe of species and total base pairs assembled in a given sample
    columns: sample, taxid, level, name, total_bp
    '''
    record[['A', 'C', 'T', 'G']] = record[['A', 'C', 'T', 'G']].astype('int')
    record['bp'] = record['A'] + record['C'] + record['T'] + record['G']
    record = record.drop(['A', 'C', 'T', 'G'], axis=1)
    taxrecord = taxrecord[['sample', 'contig', 'taxid', 'level', 'name']]
    taxrecord = taxrecord.merge(on=['sample', 'contig'], right=record)
    taxrecord = taxrecord.groupby(['sample', 'taxid', 'level', 'name'])
    taxrecord = taxrecord['bp'].sum()
    return taxrecord


def preprocess_contents():
    '''
    Process the contents large table to create the map of samples

    :input:
    =None

    :output:
    2 relating row position to samples in contents and taxonomy tables
    '''
    import xz  ## read with python-xz for speed
    import pandas as pd


    print('# working with contents:')
    fin = xz.open('contents.tsv.xz', 'rt')
    print('# getting samples mapped')
    index_sample = mapsamples(fin)
    df = pd.DataFrame.from_dict(index_sample,
                                orient='index',
                                columns=['start', 'end'])
    df = df.reset_index()
    df = df.rename({'index': 'sample'}, axis=1)
    df.to_csv('sample_map_contents.tsv.xz',
              sep='\t',
              header=True,
              index=None)
    print('# getting index dict')
    mapidx(fin,
           'sample_map_contents.tsv.xz',
           'row_indexes_contents.tsv.xz')


def preprocess_taxonomy():
    '''
    Process a large table to create the indices of samples
    in taxonomy

    :input:
    =None

    :output:
    2 relating bytes position to samples in contents and taxonomy tables
    '''
    import xz
    import pandas as pd

    print('# working with taxonomies:')
    fin = xz.open('/GMSC10/mmseqs2.lca_taxonomy.full.tsv.xz', 'rt')
    print('# getting samples mapped')
    index_sample = mapsamples(fin)
    df = pd.DataFrame.from_dict(index_sample,
                                orient='index',
                                columns=['start', 'end'])
    df = df.reset_index()
    df = df.rename({'index': 'sample'}, axis=1)
    df.to_csv('sample_map_taxonomy.tsv.xz',
              sep='\t',
              header=True,
              index=None)
    print('# getting index dict')
    mapidx(fin,
           'sample_map_taxonomy.tsv.xz',
           'row_indexes_taxonomy.tsv.xz')


def load_all():
    '''
    Load indices of samples and prepare files to 
    process by samples

    :input:
    = None

    :output:
    - posix_sample     dictionary containing row index
                       and byte position in contents file
    - posix_tax        dictionary containing row index and
                       byte position in taxonomy file
    - idx_sample       dictionary containing sample, starting and
                       ending row in contents file
    - idx_tax          dictionary containing sample, starting and
                       ending row in taxonomy file
    - fcontent         xz opened contents file
    - ftax             xz opened taxonomy file
    '''
    import xz
    import pandas as pd


    print('Load files')
    posix_sample = pd.read_table('row_indexes_contents.tsv.xz')
    posix_sample = posix_sample.set_index('sample')
    posix_tax = pd.read_table('row_indexes_taxonomy.tsv.xz')
    posix_tax = posix_tax.set_index('sample')
    fcontent = xz.open('contents.tsv.xz', 'rt')
    ftax = xz.open('/GMSC10/mmseqs2.lca_taxonomy.full.tsv.xz', 'rt')

    return [fcontent, ftax,
            posix_sample, posix_tax]


def work_merge(contents, tax, psample, ptax, ofile: str):
    '''
    Works by sample merging base-pairs per sample and contig and
    summing up if same sample and genus/species.
    
    :input:
    - contents    opened xz contents file
    - tax         opened xz taxonomy file
    - psample     dataframe with row and bytes position of contents file
    - ptax        dataframe with row and bytes position of taxonomy file
    
    :output:
    - ofile       output file containing the base-pairs by sample and genus/species    
    '''
    import pandas as pd
    print('Start preparing and merging')
    i = pd.DataFrame(columns=['sample', 'taxid', 'level', 'name', 'bp'])
    i.to_csv(ofile,
             mode='a',
             sep='\t',
             header=True,
             index=None)
    for s in psample.index:
        print(f'Process sample {s}')
        sdf1 = getdf(psample.loc[s, 'start'],
                     psample.loc[s, 'end'],
                     ['sample', 'contig', 'A', 'T', 'C', 'G'],
                     contents)
        sdf2 = getdf(ptax.loc[s, 'start'],
                     ptax.loc[s, 'end'],
                     ['sample', 'contig', 'taxid',
                      'level', 'name', 'retained',
                      'assigned', 'agreement', 'support'],
                     tax)
        sdf2 = sdf2[sdf2.level.isin(['species', 'genus'])]
        fdf = process_chunk(sdf1, sdf2)
        fdf.reset_index().to_csv(ofile,
                                 mode='a',
                                 sep='\t',
                                 header=None,
                                 index=None)


def main():
#    preprocess_contents()
#    preprocess_taxonomy()
    contents, tax, psample, ptax = load_all()
    work_merge(contents, 
               tax, 
               psample, 
               ptax, 
               './processed_contents_bp.tsv.xz')


if __name__ == '__main__':
    main()

