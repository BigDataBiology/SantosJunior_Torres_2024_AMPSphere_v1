def seqrev(seq):
    '''
    Checks if (1) sequence starts with an inverse complement of stop codon
    and if (2) ends in a stop codon.
     
    If first is true and second is true - do nothing
    If first is true and second is false - do reverse complement
    If first is false - do nothing
    '''
    if (seq[0:3] in {'TTA', 'CTA', 'TCA'}):
        if (seq[-3:] in {'TAA', 'TAG', 'TGA'}):
            pass
        else:
            return seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    return seq
    

def main():
    '''
    Plot the number of AMPs by the number of unique genes
    they are encoded by
    '''
    import lzma
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from tqdm import tqdm
    from Bio import SeqIO
    from collections import Counter


    print('Loading the genes')
    hseqs = []
    with lzma.open('data/AMPSphere_v.2022-03.fna.xz', 'rt') as infile:
        with tqdm(total=5_518_294) as pbar:
            for record in tqdm(SeqIO.parse(infile, 'fasta')):
                h = record.description.split(' | ')[1]
                s = seqrev(str(record.seq))
                hseqs.append((h, s))
                pbar.update(1)

    print('Generating a table with sequence and AMP headers')
    # AMPSphere contains 5,518,294 genes
    # testing process
    if len(hseqs) != 5_518_294:
        print('Wrong gene input process')
        exit
    
    # eliminating redundant genes
    hseqs = set(hseqs)
    df = pd.DataFrame(hseqs, columns=['AMP', 'gene'])

    print('Count the number of unique genes per amp')
    ugenes_per_amp = df.groupby('AMP').agg('size').reset_index()
    ugenes_per_amp.rename({0: '#_genes'}, axis=1, inplace=True)
    ugenes = Counter(ugenes_per_amp['#_genes'])
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
                      index=['1', '2', '>2'])
    ds = ds/1000
    
    print('Plotting unique genes')
    ds.plot.bar(legend=False, cmap='Dark2')
    plt.xlabel('Number of unique genes')
    plt.ylabel('Thousands of AMP candidates')
    plt.tight_layout()
    plt.savefig('figure_unique_genes.svg')


if __name__ == '__main__':
    main()
    
