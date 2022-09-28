def make_seqlogo(aln, oname):
    import logomaker
    from Bio import AlignIO
    from matplotlib import pyplot as plt
    # getting aligned sequences
    seqs = []
    for i in AlignIO.read(aln,
                          'fasta'):
        seqs.append(str(i.seq))
    # convert to matrix
    df = logomaker.alignment_to_matrix(seqs)
    # transform
    df = logomaker.transform_matrix(df,
                                    from_type='counts',
                                    to_type='information')
    # make logo
    seqlogo = logomaker.Logo(df,
                            font_name='DejaVu Sans',
                            color_scheme='skylign_protein',
                            vpad=.1,
                            width=.8)
    # style using Logo methods
    seqlogo.style_xticks(anchor=0,
                         spacing=5,
                         rotation=45)
    # style using Axes methods
    seqlogo.ax.set_ylabel('information (bits)')
    seqlogo.ax.set_xlim([-1,
                         len(df)])
    plt.tight_layout()
    plt.savefig(oname)
    plt.close()


def batch_seqlogo():
    import os
    from glob import glob
    data_folder = 'analysis/families/aln'
    analysis_folder = 'analysis/families/seqlogo'
    os.makedirs(analysis_folder, exist_ok=True)
    for f in glob(f'{data_folder}/*.aln'):
        infile = f.split('/')[-1]
        ofile = infile.replace('.aln', '.svg')
        make_seqlogo(f'{data_folder}/{infile}',
                     f'{analysis_folder}/{ofile}')
        print(ofile)             

