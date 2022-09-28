def fix_multigenes():
    from glob import glob
    from Bio import SeqIO
    
    for sample in glob('analysis/*.fa'):
        seqs = []
        for record in SeqIO.parse(sample, 'fasta'):
            seqs.append((record.id, str(record.seq)))
        if len(seqs) > 1:
            L = [len(x) for _, x in seqs]
            n = max(L)
            for l, (h, s) in zip(L, seqs):
                if l == n:
                    with open(sample, 'wt') as ofile:
                        ofile.write(f'>{h}\n{s}\n')


def makeone():
    from glob import glob
    from Bio import SeqIO
    
    seqs = []
    for sample in glob('analysis/*.fa'):
        s = sample.replace('.fa', '')    
        for record in SeqIO.parse(sample, 'fasta'):
            seqs.append((s, str(record.seq)))
    
    with open('analysis/p1_genes.fasta', 'wt') as ofile:
        for h, s in seqs:
            ofile.write(f'>{h}\n{s}\n')


def main():
    fix_multigenes()
    makeone()


if __name__ == '__main__':
    main()
    
