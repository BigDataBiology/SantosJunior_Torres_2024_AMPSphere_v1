def homologs_search_per_db():
    import pandas as pd
    
    infolder = 'analysis/homologs'
    
    dramp = pd.read_table(f'{infolder}/dramp_candidates.txt', header=None)
    gmgc = pd.read_table(f'{infolder}/gmgc_candidates.txt', header=None)
    smprot = pd.read_table(f'{infolder}/SmProt_candidates.txt', header=None)
    starpep = pd.read_table(f'{infolder}/starPepDB_candidates.txt', header=None)
    STsorfs = pd.read_table(f'{infolder}/STsORFs_candidates.txt', header=None)
    
    def formatname(x):
        x = str(x).zfill(6)
        return '_'.join([x[idx:idx + 3]
                         for idx, val in enumerate(x) if idx % 3 == 0])
    
    all_amps = [ f'AMP10.{formatname(x)}' for x in range(863499)]
    df = pd.DataFrame(all_amps, columns=['AMPSphere'])
    df[['DRAMP', 'GMGCv1', 'SmProtv2', 'StarPepDB45k', 'STsORFs']] = False
    
    df = df.set_index('AMPSphere')
    df.loc[dramp[0], 'DRAMP'] = True
    df.loc[gmgc[0], 'GMGCv1'] = True
    df.loc[smprot[0], 'SmProtv2'] = True
    df.loc[starpep[0], 'StarPepDB45k'] = True
    df.loc[STsorfs[0], 'STsORFs'] = True
    df = df.reset_index()
    
    df.to_csv(f'{infolder}/homologs_table.tsv',
              sep='\t',
              header=True,
              index=None)
    
    df.DRAMP = ['DRAMP' if x == True else '' for x in df.DRAMP]
    df.GMGCv1 = ['GMGCv1' if x == True else '' for x in df.GMGCv1]
    df.SmProtv2 = ['SmProtv2' if x == True else '' for x in df.SmProtv2]
    df.StarPepDB45k = ['StarPepDB45k' if x == True else '' for x in df.StarPepDB45k]
    df.STsORFs = ['STsORFs' if x == True else '' for x in df.STsORFs]
    
    d = df.set_index('AMPSphere').apply(list, axis=1)
    d = d.reset_index()
    d = d.rename({0: 'homologs'}, axis=1)
    
    d.homologs = d['homologs'].apply(lambda x: set(x))
    d.homologs = d.homologs.apply(lambda x: ', '.join(x))
    
    d.to_csv(f'{infolder}/homologs_list.tsv',
             sep='\t',
             header=True,
             index=None)


def enrichment_search():
    import pandas as pd
    
    hq = pd.read_table('data/high_quality_candidates.txt', header=None)
    q = pd.read_table('data/quality_candidates.txt', header=None)
    qc = pd.concat([hq, q])
    
    homologs = pd.read_table('analysis/homologs/homologs_table.tsv.txt')
    homologs = homologs.set_index('AMPSphere')
    db_totals = homologs.sum(axis=0)
    
    dbs = ['DRAMP', 'GMGCv1', 'SmProtv2', 'StarPepDB45k', 'STsORFs']
    
    ratio_ampsphere = len(qc)*100/len(homologs)
    
    enrichment_results = []
    for i in dbs:
        r = homologs.loc[qc[0], i].sum(axis=0)
        new = [i, r]
        r = r*100 / db_totals[i]
        new.append(db_totals[i])
        new.append(r)
        new.append(r/ratio)
        new.append(r/ratio_ampsphere)
        enrichment_results.append(new)
    
    homologs_assessment = pd.DataFrame(enrichment_results,
                                       columns=['database', 'hq_peptides',
                                                'total_homologs', 'ratio',
                                                'enrichment_from_AMPSphere'])
    
    homologs_assessment.to_csv('analysis/homologs/homologs_assessment_per_db.tsv',
                               sep='\t',
                               header=True,
                               index=None)


def plot_qual_homologs():
    import pandas as pd
    import seaborn as sns
    from matplotlib import pyplot as plt
    
    qualassess = pd.read_table('data/quality_assessment.tsv')
    homologs = pd.read_table('analysis/homologs/homologs_table.tsv')
    homologs.rename({'AMPSphere': 'AMP'}, axis=1, inplace=True)
    
    df = qualassess.merge(on='AMP', right=homologs)
    df = df.replace(['Passed', 'Failed', 'Not tested'], [True, False, False])
    
    tests = list(qualassess.columns)[1:]
    databases = list(homologs.columns)[1:]
    
    newdf = pd.DataFrame()
    for i in databases:
        n = df[df[i] == True].sum(axis=0).loc[tests]
        d = df[df[i] == True].sum(axis=0).loc[i]
        newdf[i] = n * 100 / d
    
    newdf = newdf.astype('float16')
    newdf.index = ['Antifam', 'RNAcode',
                   'Metaproteomes',
                   'Metatranscriptomes',
                   'Terminal position']
    newdf.to_csv('analysis/homologs/quality_by_db.tsv',
                 sep='\t',
                 header=True,
                 index=True)  
    
    sns.clustermap(newdf, cmap='YlOrBr')
    plt.savefig('analysis/homologs/quality_over_homologs.svg')
    plt.close()
    
    
def dha():
    homologs_search_per_db()
    enrichment_search()
    plot_qual_homologs()


