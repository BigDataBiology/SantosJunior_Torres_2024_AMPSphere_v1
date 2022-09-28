def create_consensus():
    from glob import glob
    from Bio import AlignIO, SeqIO
    from Bio.Align import AlignInfo
    from Bio.SeqRecord import SeqRecord

    seqs = []
    for i in glob('analysis/families/aln/*.aln'):
        aln = AlignIO.read(i, 'fasta')
        saln = AlignInfo.SummaryInfo(aln)
        fam = i.split('/')[-1].replace('.aln', '')
        record = SeqRecord(saln.dumb_consensus(0.75),
                           id=fam,
                           name=fam,
                           description='consensus at 75% of conservation per site')
        seqs.append(record)

    with open('analysis/consensus_seqs.fasta', 'w') as output_handle:
        SeqIO.write(seqs, output_handle, 'fasta')


def align_consensus():
    from Bio.Align.Applications import MafftCommandline
    mafft_exe = '/usr/bin/mafft'
    in_file = 'analysis/consensus_seqs.fasta'
    mafft_cline = MafftCommandline(mafft_exe, input=in_file)

    stdout, stderr = mafft_cline()
    with open('analysis/aligned_consensus.fasta', 'w') as handle:
        handle.write(stdout)


def trees():
    from .phylogen import treebuilder, draw_ascii
    treebuilder('analysis/aligned_consensus.fasta',
                'analysis/aligned_consensus.nwk',
                'WAG+CAT',
                1000)

    draw_ascii('analysis/aligned_consensus.nwk',
               'analysis/aligned_consensus.ascII',
               column_width=80)


def organizer():
    import os
    infiles = ['aligned_consensus.fasta',
               'aligned_consensus.nwk',
               'aligned_consensus.ascII',
               'consensus_seqs.fasta']
    ofolder = 'analysis/families'
    ifolder = 'analysis'
    for i in infiles:
        os.rename(f'{ifolder}/{i}',
                  f'{ofolder}/{i}')


def canalyzer():
    create_consensus()
    align_consensus()
    tress()
    organizer()
    
