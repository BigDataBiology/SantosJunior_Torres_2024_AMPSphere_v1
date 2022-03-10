def genes():
    '''
    Retrieve genes from AMPSphere
    and generate an output as fna file
    '''
    import lzma
    import pandas as pd
    from Bio import SeqIO
    
    data_folder = 'data/'
    analysis_folder = 'analysis/'
    
    corrlist = f'{analysis_folder}/AMPsphere_GMSC_correspondence.tsv.gz' 
    corrlist = pd.read_table(corrlist, sep='\t', header='infer')
    
    # generating correspondence list of access and genes
    corrlist = corrlist[['accession', 'genes']]
    corrlist = corrlist.set_index('accession')
    corrlist['genes'] = [x.split(',') for x in corrlist.genes]

    dictresource = dict()
    for f in corrlist.itertuples():
        for x in f[1]: dictresource[x] = f[0]

    # generating output
    fout = f'{analysis_folder}/AMPSphere_v.2021-03.fna.xz'
    
    fin = f'{data_folder}/gmsc_genes.fna.xz'
    fin = lzma.open(fin, 'rt', encoding='utf-8')
    
    with lzma.open(fout, 'wt', encoding='utf-8') as dbout:
        for record in SeqIO.parse(fin,
                                  'fasta'):
            header = f'>{record.id} | {dictresource[record.id]}\n'
            print(header.strip())
            dbout.write(f'{header}{str(record.seq)}\n')

    fin.close()
    
