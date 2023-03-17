def ugenes_plot():
    '''
    Plot the number of AMPs by the number of unique genes
    they are encoded by
    '''
    import lzma
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    
    from Bio import SeqIO
    from collections import Counter


    print('Loading the genes')
    headers, seqs = [], []
    with lzma.open('data/AMPSphere_v.2022-03.fna.xz', 'rt') as infile:
        for record in SeqIO.parse(infile, 'fasta'):
            headers.append(record.description.split(' | ')[1])
            seqs.append(str(record.seq))

    print('Generating a table with sequence and AMP headers')
    df = pd.DataFrame(np.array([headers, seqs]).T, columns=['AMP', 'gene'])
    # eliminating redundant genes
    df = df.drop_duplicates()

    print('Count the number of unique genes per amp')
    ugenes_per_amp = dict(Counter(df.AMP))
    ugenes = Counter(ugenes_per_amp.values())
    ds = pd.DataFrame.from_dict(ugenes,
                                orient='index',
                                columns=['number_of_amps'])
    ds = ds.reset_index()
    ds = ds.rename({'index': 'number_of_ugenes'},
                   axis=1)

    # sorting by number of ugenes
    ds = ds.sort_values(by='number_of_ugenes')
    
    # counting number of amps with 3 ugenes or more
    amps = ds[ds.number_of_ugenes >= 3]['number_of_amps'].sum()

    # organizing dataframe
    ds = ds.set_index('number_of_ugenes')
    ds = pd.DataFrame([ds.loc[1]['number_of_amps'],
                       ds.loc[2]['number_of_amps'],
                       amps],
                      index=['1', '2', '>3'])

    print('Plotting unique genes')
    ds.plot.bar(legend=False, cmap='Dark2')
    plt.xlabel('Number of unique genes')
    plt.ylabel('AMP candidates')
    plt.tight_layout()
    plt.savefig('figure_S1b_unique_genes.svg')

