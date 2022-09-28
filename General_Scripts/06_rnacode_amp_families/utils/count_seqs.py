def format_amp(ampgene):
    ampgene = ampgene.split('_')
    ampgene = ampgene[:-1]
    ampgene = '_'.join(ampgene)
    return ampgene


def lenaln(f):
    from Bio import AlignIO
    a = AlignIO.read(f, 'clustal')
    ids = [format_amp(j.id) for j in a]
    ids = set(ids)
    l = len(a)
    return [','.join(ids), l]
   
    
def measure():
    from os.path import exists
    from glob import glob
    records = []
    for f in glob('*.aln'):
        print(f)
        if exists(f):
            mlist = lenaln(f)
            f = f.replace('.aln', '')
            mlist.append(f)
            records.append(mlist)
    return records
    

def dfparse(records):
    import pandas as pd
    data = pd.DataFrame(records,
                        columns=['amps',
                                 'genes',
                                 'cluster'])
    data = data.sort_values(by=['genes', 'cluster'])
    data.to_csv('analysis/file_order.tsv',
                sep='\t',
                header=True,
                index=None)
    return data


def count_seqs():
    records = measure()
    data = dfparse(records)
    return data
    
