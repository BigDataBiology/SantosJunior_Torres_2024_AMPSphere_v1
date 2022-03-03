def getnumber(header):
    '''
    Extract number from header string
    '''
    header = header.split('GMSC10.SMORF.')[1]
    header = header.replace('_', '')

    return int(header)


def format_progenomes(ampdict: dict, line: str):
    '''
    Takes the gene coordinates in the genome
    and return a set of variables
    '''
    line = line.strip().split(':')
    line[1] = line[1].split('-')  # start-stop
    start = line[1][0]
    stop = line[1][1]
    contig = line[0].split('_')[0]
    genome = '.'.join(line[0].split('.')[0:2])
    gene = genes[i]
    amp = ampdict[gene]

    return (amp, gene, genome, contig, start, stop)


def ampsphere2progenomes():
    '''
    Associate AMPSphere to progenomes2 database
    '''
    import gzip
    import lzma
    import pandas as pd
    from preprocess import create_GMGCgenes
    from preprocess import genes2amp

    # select only genes from progenomes
    genes = create_GMGCgenes()
    for i in genes:
        n = getnumber(i)
        if n > 34_617_405:  # metagenomic gmsc genes start
            genes.remove(i)
            
    # available for the Zenodo repository
    ampdict = genes2amp()

    out = gzip.open('data/AMPSphere_ProGenomes2.tsv.gz',
                    'wt',
                    encoding='utf-8')

    out.write('AMP\tGMSC10\tgenome\tcontig\tstart\tstop\n')

    # this table is internal
    i = 0
    with gzip.open('data/GMSC10.ProGenomes2.coords.txt.gz',
                   'rt',
                   encoding='utf-8') as d:
        for j, line in enumerate(d):
            if i < len(genes) and j < 34_617_405:
                n = getnumber(genes[i])
                if j == int(n):
                    (amp, gene, genome, contig, start, stop) = format_progenomes(ampdict, line)
                    out.write(f'{amp}\t{gene}\t{genome}\t{contig}\t{start}\t{stop}\n')
                    i += 1

    out.close()
    
