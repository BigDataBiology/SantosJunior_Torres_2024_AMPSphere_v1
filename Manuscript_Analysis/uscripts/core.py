import gzip
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter


def classificationprop(x):
    '''
    Classify the prevalence percent
    into different categories
    '''
    if x >= 90: return 'core'
    if x >= 50: return 'shell'
    if x < 50: return 'accessory'


def progenomes_prime():
    '''
    Determines the size of specI clusters
    '''
    df = pd.read_table('data/pgenomes_samples.tsv.gz',
                       sep='\t',
                       header=None,
                       names=['specI', 'genome'])

    spec_size = dict(Counter(df['specI']))

    return spec_size


def linksource():
    '''
    Loads the sources table and return 
    a smaller version of itself only with meaningful
    columns to this analysis
    '''
    # linking to the source resource
    refdata = pd.read_table('data/gmsc_amp_genes_envohr_source.tsv.gz',
                            sep='\t',
                            header='infer')

    refdata = refdata[refdata.is_metagenomic == False]
    refdata = refdata[~(refdata.specI.isna())]
    refdata = refdata[['amp','sample','specI']]
    refdata = refdata.sort_values(by=['amp', 'sample'])
    refdata = refdata.drop_duplicates()

    return refdata


def famp_analysis(spec_size, refdata):
    '''
    Get the families and classify their occurrence into
    core, shell or strain-specific (accessory)
    '''
    
    # creating a dictionary of AMPs to families
    spheres = pd.read_table('data/SPHERE_v.2022-03.levels_assessment.tsv.gz',
                            sep='\t',
                            header='infer')
    spheres = spheres[['AMP accession', 'SPHERE_fam level III']]
    spheres = spheres.rename({'AMP accession': 'amp',
                              'SPHERE_fam level III': 'families'},
                             axis=1)
    
    refdata = refdata.merge(on='amp', right=spheres)
    
    ## analyzing families
    fams = refdata[['families', 'sample', 'specI']]
    fams = fams.drop_duplicates()
    fams = fams.drop('sample', axis=1)
    fams = fams.sort_values(by=['families', 'specI'])

    # counting clusters per family
    ofile = gzip.open('families_all.count_core.tsv.gz', 'wt', encoding='utf-8')
    ofile.write('family\tspecI\tcounts\ttotal\tproportion\tclassification\n')

    for df in fams.groupby('families'):
        f = df[0]
        D = dict(Counter(df[1].specI))
        for k, v in D.items():
            if spec_size.get(k, 0) >= 10:
                nv = v * 100 / spec_size[k]
                total = spec_size[k]
                ofile.write(f'{f}\t{k}\t{v}\t{total}\t{nv}\t{classificationprop(nv)}\n')

    ofile.close()

    ## analyzing amps
    amps = refdata[['amp', 'sample', 'specI']]
    amps = amps.drop_duplicates()
    amps = amps.drop('sample', axis=1)

    # counting clusters per family
    ofile = gzip.open('amps_all.count_core.tsv.gz', 'wt', encoding='utf-8')
    ofile.write('amp\tspecI\tcounts\ttotal\tproportion\tclassification\n')

    for df in amps.groupby('amp'):
        f = df[0]
        D = dict(Counter(df[1].specI))
        for k, v in D.items():
            if spec_size.get(k, 0) >= 10:
                nv = v * 100 / spec_size[k]
                total = spec_size[k]
                ofile.write(f'{f}\t{k}\t{v}\t{total}\t{nv}\t{classificationprop(nv)}\n')

    ofile.close()


def core():
    print('Get size of specI clusters')
    spec_size = progenomes_prime()
    print('Link sources')
    refdata = linksource()
    print('Analyzing families')
    famp_analysis(spec_size, refdata)
    
    
def main():
    '''
    Select quality families and AMPs
    and plot the bar chart with them
    including the info about the core,
    shell and accessory classifications
    '''

    # perform core analysis
    core()
    
    # load data
    data = pd.read_table('data/SPHERE_v.2022-03.levels_assessment.tsv.gz')
    
    # define minimum size of a family as 8 AMPs to be plotted 
    minsizedfam = [k for k,v in Counter(data['SPHERE_fam level III']).items() if v >= 8]

    # filter families by size and quality
    data = pd.read_table('data/quality_families.txt.gz')
    qualfam = data[data.total >= 8]
    qualfam = qualfam[(qualfam.experimental_evidence == True) | (qualfam.perc >= 75.0)]
    qualfam = qualfam.family.tolist()

    # getting quality-controlled AMPs
    data = pd.read_table('data/quality_candidates.txt.gz', header=None)
    data2 = pd.read_table('data/high_quality_candidates.txt.gz', header=None)
    qualamp = pd.concat([data, data2])[0].tolist()

    # loading core analysis for families
    d = pd.read_table('families_all.count_core.tsv.gz', sep='\t', header='infer')

    # creating set for all families with minimum size
    dfam = d[d.family.isin(minsizedfam)]
    dfam = dfam.classification
    dfam = Counter(dfam)    
    dfam = pd.DataFrame.from_dict(dfam,
                                  orient='index',
                                  columns=['All_families']).T

    # creating set for families with minimum size and quality
    dqcfam = d[d.family.isin(qualfam)]
    dqcfam = dqcfam.classification
    dqcfam = Counter(dqcfam)
    dqcfam = pd.DataFrame.from_dict(dqcfam,
                                    orient='index',
                                    columns=['HQ_families']).T


    # loading core analysis for families
    d = pd.read_table('amps_all.count_core.tsv.gz', sep='\t', header='infer')

    # creating set for all AMPs
    dfamp = Counter(d.classification)
    dfamp = pd.DataFrame.from_dict(dfamp,
                                   orient='index',
                                   columns=['All_AMPs']).T

    # creating set for AMPs with quality
    dqcamp = d[d.amp.isin(qualamp)]
    dqcamp = dqcamp.classification
    dqcamp = Counter(dqcamp)    
    dqcamp = pd.DataFrame.from_dict(dqcamp,
                                    orient='index',
                                    columns=['HQ_AMPs']).T

    # concatenating separated dataframes
    df = pd.concat([dfamp, dqcamp, dfam, dqcfam]).T

    # calculating percent to normalize the measures
    dnew = df * 100 / df.sum()

    # plotting data
    dnew.T.plot.barh(stacked=True, cmap='Dark2')
    plt.xlim(50,100)
    plt.xlabel('Percent (%)')
    plt.ylabel('Classification')
    plt.tight_layout()
    plt.savefig('core_analysis.svg')
    plt.close()


if __name__ == '__main__':
    main()
    
