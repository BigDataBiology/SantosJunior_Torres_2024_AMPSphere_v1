def openfasta(infile):
    from Bio import SeqIO
    samplename = infile.split('.')[0]
    samplename = samplename.split('/')[-1]
    if infile.endswith('.gz'):
        import gzip
        infile = gzip.open(infile,
                           'rt')       
    hs = []
    for idx, record in enumerate(SeqIO.parse(infile, 'fasta')):
        header = str(idx).zfill(12)
        header = samplename+'_'+header
        hs.append((header, str(record.seq)))
    return hs


