def helical_wheel():
    '''
    Generate helical wheels
    for each sequence in
    a fasta file
    '''
    import sys
    import os
    import gzip
    from Bio import SeqIO
    from modlamp.plot import helical_wheel
    import matplotlib.pyplot as plt

    analysis_folder='analysis/'

    hf = f'{analysis_folder}/helical_wheels/'
    os.makedirs(hf,
                exist_ok=True)

    infile = gzip.open(f'{analysis_folder}/AMPSphere_v.2021-03.faa.gz',
                       'rt',
                       encoding='utf-8')

    for index, record in enumerate(SeqIO.parse(infile,"fasta")):
            file_name = f'{hf}/helicalwheel_{str(record.id)}.svg'
            helical_wheel(str(record.seq),
                          moment=True,
                          colorcoding='amphipathic',
                          lineweights=True,
                          filename=file_name,
                          seq=False)
            plt.close()
            print(file_name)

    infile.close()

