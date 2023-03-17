def rnacode(in_file, out_file):
    import subprocess
    subprocess.call(['RNAcode',
                     '-o',
                     out_file,
                     '-t',
                     '--best-only',
                     '--stop-early',
                     '-n', '1000',
                     '-p', '0.05',
                     in_file])


def batch_process():
    import os
    from .count_seqs import count_seqs
    os.makedirs('analysis/alignments', exist_ok=True)
    os.makedirs('analysis/rnacode_out', exist_ok=True)
    data = count_seqs()
    for in_file in data.cluster:
        rnacode(f'{in_file}.aln',
                f'analysis/rnacode_out/{in_file}.rnacode')
        os.rename(f'{in_file}.aln',
                  f'analysis/alignments/{in_file}.aln')

