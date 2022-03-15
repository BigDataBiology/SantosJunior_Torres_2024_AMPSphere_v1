def amplen():
    '''
    AMP features distribution
    in AMPSphere
    '''
    import gzip
    import pandas as pd
    import matplotlib.pyplot as plt

    from Bio import SeqIO
    
    
    print('... loading features')
    data = pd.DataFrame()
    Lamp = []
    with gzip.open('data/AMPSphere_v.2022-03.faa.gz', 'rt', encoding='utf-8') as ifile:
        for record in SeqIO.parse(ifile, 'fasta'):
            Lamp.append(len(record.seq))
            
    data['length'] = Lamp

    # preparing the histogram of isolectric points
    data.length.hist(bins=100, grid=False)
    plt.xlim(5,100)
    plt.ylabel('Number of AMP candidates')
    plt.xlabel('Length in residues')
    plt.savefig('figure_S1c_histogram_AMPs_length.svg')
    plt.close()

