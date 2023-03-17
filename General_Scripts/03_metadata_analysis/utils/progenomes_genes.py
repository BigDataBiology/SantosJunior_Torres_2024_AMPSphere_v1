def getnumber(header):
    '''
    Extract number from header string
    '''
    header = header.split('GMSC10.SMORF.')[1]
    header = header.replace('_', '')

    return int(header)


def format_progenomes(line: str):
    '''
    Takes the gene coordinates in the genome
    and return a set of variables
    '''
    line = line.strip().split(':')
    line[1] = line[1].split('-')  # start-stop
    start = line[1][0]
    stop = line[1][1]
    contig = line[0].split('.')[2]
    genome = '.'.join(line[0].split('.')[0:2])

    return (genome, contig, start, stop)


def create_GMGCgenes():
    '''
    Creates a list of GMSC genes
    '''
    import pandas as pd
    
    data = pd.read_table('data/AMPSphere_v.2022-03.origin_samples.tsv.gz', sep='\t', header='infer')

    gmscgenes = []
    for i in data['genes']: gmscgenes += i.split(',')    

    gmscgenes.sort()
        
    ampdict = dict()
    for i in data[['accession', 'genes']].itertuples():
        for j in i[2].split(','): ampdict[j] = i[1]
    
    return gmscgenes, ampdict


def ampsphere2progenomes():
    '''
    Associate AMPSphere to progenomes2 database
    '''
    import gzip
    import lzma
    import pandas as pd
    from preprocess import genes2amp

    # select only genes from progenomes
    genes, ampdict = create_GMGCgenes()

    # metagenomic gmsc genes start
    genes = dict()
    for k in gmscgenes:
        if getnumber(k) <= 34_617_405:
            genes[getnumber(k)] = k
    
    out = gzip.open('data/AMPSphere_ProGenomes2.tsv.gz',
                    'wt',
                    encoding='utf-8')

    out.write('AMP\tGMSC10\tgenome\tcontig\tstart\tstop\n')

    # this table is internal
    with gzip.open('data/GMSC10.ProGenomes2.coords.txt.gz',
                   'rt',
                   encoding='utf-8') as d:
        for j, line in enumerate(d):
            if j in genes:
                (genome, contig, start, stop) = format_progenomes(line)
                gene = genes[j]
                amp = ampdict[gene]
                out.write(f'{amp}\t{gene}\t{genome}\t{contig}\t{start}\t{stop}\n')

    out.close()
    
