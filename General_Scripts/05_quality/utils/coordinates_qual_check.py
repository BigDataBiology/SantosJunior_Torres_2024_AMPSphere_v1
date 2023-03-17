def returngenes():
    import lzma
    from Bio import SeqIO
    seqs = dict()
    for record in SeqIO.parse(lzma.open('data/AMPSphere_v.2022-03.fna.xz',
                                        'rt'),
                              'fasta'):
        x = record.description.split(' | ')
        seqs[x[0]] = x[1]
    return seqs


def recovergmsc(x, seqs):
    for k, v in seqs.items():
        if v == x: yield k


def sampling():
    import pandas as pd
    data = pd.read_table('analysis/quality_assessment.tsv')
    data = data[['AMP', 'Coordinates']]
    data = data[data.Coordinates == 'Passed']['AMP']
    s1 = data.sample(100).tolist()
    s2 = data[(~data.isin(s1))].sample(100).tolist()
    s3 = data[(~data.isin(s1+s2))]
    return s1, s2, s3


def rcoordinates(samples):
    import pandas as pd
    df = pd.DataFrame()
    for chunk in pd.read_table('data/coordinates_test_passed.tsv.gz',
                               sep='\t',
                               header=None,
                               names=['gene', 'status'],
                               chunksize=25_000_000):
        chunk = chunk[chunk.gene.isin(samples)]
        df = pd.concat([df, chunk])
        print(len(df))
    return df
    

def makeassoc(s1, s2, s3):
    assoc = dict()
    for i in s1: assoc[i] = 'rep1'
    for i in s2: assoc[i] = 'rep2'
    for i in s3: assoc[i] = 'rep3'
    return assoc
    
    
def get_status(df, s1, s2, s3, seqs):
    import pandas as pd
    from collections import Counter
    df['amp'] = [seqs[x] for x in df.gene]
    table = pd.pivot_table(df, index=['amp'], columns=['status'], aggfunc='count')
    table = table.fillna(0)
    table.columns = ['F', 'T']
    table['%T'] = table['T']*100/(table['T'] + table['F'])
    names = makeassoc(s1, s2, s3)
    table['replicate'] = [names[x] for x in table.index]
    return table.reset_index()
    

def main():
    from itertools import chain
    seqs = returngenes()
    s1, s2, s3 = sampling()
    sample1 = [list(recovergmsc(x, seqs)) for x in s1]
    sample2 = [list(recovergmsc(x, seqs)) for x in s2]
    sample3 = [list(recovergmsc(x, seqs)) for x in s3]
    sample1 = list(chain.from_iterable(sample1))
    sample2 = list(chain.from_iterable(sample2))
    sample3 = list(chain.from_iterable(sample3))
    samples = sample1 + sample2 + sample3
    df = rcoordinates(samples)
    table = get_status(df, s1, s2, s3, seqs)
    table.to_csv('analysis/quality_coordinates_check.tsv', sep='\t', header=True, index=None)
    ndf = table.drop(['amp', '%T'], axis=1).groupby('replicate').agg('sum')
    ndf['T_pct'] = ndf['T'] * 100 / (ndf.F + ndf['T'])
    ndf.to_csv('analysis/quality_coordinates_summary_check.tsv', sep='\t', header=True, index=True)
    

