def motif():
    motif = ['tle target', 'ion-pair-pi', 'vinculin',
             'gly-central–symmetrical', 'leulystrp',
             'cbp-a', 'dock_MAPK', 'guanine nucleotide exchange',
             'kld12', 'γ-core motif', 'golgi_lock',
             'arginine finger motif', 'cargo', 'walker',
             'bact7', 'cspa', 'kwarn', 'lfcinb', 'atcun',
             'pest-motif', 'siderophore', 'eno',
             'autoproteolytic cleavage motif', 'chemotaxis', 'spsa',
             'ngr motif', 'gly zipper', 'n-myristoylation/s-palmitoylation motif',
             'yada', 'ABC_motif', 'adherence', 'capsule', 'rgd',
             'lipid binding motif', 'gala', 'mtpm', 'rana box', 'type_iv',
             'bfii', 'bact5', 'caax box']
    return motif
    
    
def load_annotations():
    import pandas as pd
    data = pd.read_table('analysis/AMPSphere_v.2022-03.annotation.tsv',
                         sep='\t',
                         header='infer')
    data = data.fillna('n.a.')
    data['motif_match'] = [x.split('|') for x in data.motif_match]
    data['AA_absent'] = [x.split('|') for x in data.AA_absent]
    data['AA_rich'] = [x.split()[0] for x in data.AA_rich]
    data['log2(gene_prob)'] = [2**x for x in data['log2(gene_prob)']]
    data['log2(gene_prob)'] = 863498 * data['log2(gene_prob)']
    data.rename({'log2(gene_prob)': 'prob_at_random'},
                axis=1,
                inplace=True) 
    return data
    
    
def load_clusters():
    import pandas as pd
    data = pd.read_table('data/SPHERE_v.2022-03.levels_assessment.tsv.gz',
                         sep='\t', 
                         header='infer')
    data = data[['AMP accession', 'SPHERE_fam level III']]
    data = data.rename({'AMP accession': 'id',
                        'SPHERE_fam level III': 'family'},
                       axis=1)
    return data


def merge():
    anno = load_annotations()
    fam = load_clusters()
    return anno.merge(on='id', right=fam)


def famsize():
    import pandas as pd
    from collections import Counter
    fam = load_clusters()
    fam = pd.DataFrame.from_dict(Counter(fam['family']),
                                 orient='index',
                                 columns=['fam_size'])
    fam = fam.reset_index()
    fam = fam.rename({'index': 'family'},
                     axis=1)
    fam = fam.sort_values(by='family')
    return fam
    

def sumlists(d):
    from collections import Counter
    neo = []
    for i in d: neo = neo + i
    return dict(Counter(neo))


def df_generator(df):
    import pandas as pd
    from collections import Counter
    fam = famsize()
    fam = fam[fam.fam_size == 1]['family']

    neo = [pd.DataFrame([],
                       columns=motif())]
    naa = [pd.DataFrame([], 
                       columns=['A', 'C', 'D', 'E',
                                'F', 'G', 'H', 'I',
                                'K', 'L', 'M', 'N',
                                'P', 'Q', 'R', 'S',
                                'T', 'V', 'W', 'Y'])]
    nar = [pd.DataFrame([],
                       columns=['A', 'C', 'D', 'E',
                                'F', 'G', 'H', 'I',
                                'K', 'L', 'M', 'N',
                                'P', 'Q', 'R', 'S',
                                'T', 'V', 'W', 'Y'])]

    for record in df.groupby('family'):
        print(record[0])
        neo.append(pd.DataFrame.from_dict(sumlists(record[1]['motif_match']),
                                          orient='index',
                                          columns=[record[0]]).T)
        naa.append(pd.DataFrame.from_dict(sumlists(record[1]['AA_absent']),
                                          orient='index',
                                          columns=[record[0]]).T)
        nar.append(pd.DataFrame.from_dict(Counter(record[1]['AA_rich']),
                                          orient='index',
                                          columns=[record[0]]).T)
    print('Concat the enriched amino acids')
    nar = pd.concat(nar)
    print('Concat the motifs')
    neo = pd.concat(neo)
    print('Concat the absent amino acids')
    naa = pd.concat(naa)
    return neo, naa, nar


def plot_bar(df, figtitle, ylab, xlab):
    import matplotlib.pyplot as plt
    df = df.set_index('family')
    df = df.astype('int')
    df = df.sum(axis=0)
    df = df.sort_values()
    df.plot.bar()
    plt.ylabel(ylab)
    plt.xlabel(xlab)
    plt.tight_layout()
    plt.savefig(f'analysis/{figtitle}.svg')
    plt.close()
    

def fixdfs(neo, naa, nar):
    neo = neo.fillna('0')
    neo = neo.reset_index()
    neo = neo.rename({'index': 'family'}, axis=1)
    neo.to_csv('analysis/motif_seqs_per_family.tsv',
               sep='\t',
               header=True,
               index=None)    
    naa = naa.fillna('0')
    naa = naa.reset_index()
    naa = naa.rename({'index': 'family'}, axis=1)
    naa.to_csv('analysis/absent_aa_seqs_per_family.tsv',
               sep='\t',
               header=True,
               index=None)    
    nar = nar.fillna('0')
    nar = nar.reset_index()
    nar = nar.rename({'index': 'family'}, axis=1)
    nar.to_csv('analysis/rich_aa_seqs_per_family.tsv',
               sep='\t',
               header=True,
               index=None)
    

def plot_dfs():
    import pandas as pd
    neo = pd.read_table('analysis/motif_seqs_per_family.tsv')
    naa = pd.read_table('analysis/absent_aa_seqs_per_family.tsv')
    nar = pd.read_table('analysis/rich_aa_seqs_per_family.tsv')
    plot_bar(neo,
             'motif_seqs',
             'AMPs',
             'Motifs')
    plot_bar(nar,
             'enrich_aa_seqs',
             'AMPs',
             'Enriched amino acids')
    plot_bar(naa,
             'absent_aa_seqs',
             'AMPs',
             'Absent amino acids')


def get_motifs(row):
    row = row[1]
    row = row[row > 75]
    row = row.index.tolist()
    return row
    

def normdfs():        
    import pandas as pd
    fam = famsize()
    neo = pd.read_table('analysis/motif_seqs_per_family.tsv')
    naa = pd.read_table('analysis/absent_aa_seqs_per_family.tsv')
    nar = pd.read_table('analysis/rich_aa_seqs_per_family.tsv')
    neo = neo.merge(on='family', right=fam)
    neo = neo.set_index('family')
    neo = neo.iloc[:,:-2].div(neo['fam_size'], axis=0)*100   
    naa = naa.merge(on='family', right=fam)
    naa = naa.set_index('family')
    naa = naa.iloc[:,:-2].div(naa['fam_size'], axis=0)*100
    nar = nar.merge(on='family', right=fam)
    nar = nar.set_index('family')
    nar = nar.iloc[:,:-2].div(nar['fam_size'], axis=0)*100
    neo.to_csv('analysis/motifs_norm_family.tsv',
               sep='\t',
               header=True,
               index=True)
    naa.to_csv('analysis/absent_aa_norm_family.tsv',
               sep='\t',
               header=True,
               index=True)
    nar.to_csv('analysis/rich_aa_norm_family.tsv',
               sep='\t',
               header=True,
               index=True)


def countfams():
    import pandas as pd
    import matplotlib.pyplot as plt
    qual = pd.read_table('data/quality_families.txt')
    qual = qual[qual.total >= 8]
    qual = qual[(qual.experimental_evidence == True) | (qual.perc >= 75)]
    qual = qual['family'].tolist()
    print('Processing motifs')
    neo = pd.read_table('analysis/motifs_norm_family.tsv')
    neo = neo.set_index('family')
    mlist = []
    for row in neo.iterrows():
        n = get_motifs(row)
        if len(n) > 0:
            mlist.append([row[0], n])
            
    mlist = pd.DataFrame(mlist,
                         columns=['family',
                                  'motifs'])
    mlist = mlist[mlist.family.isin(qual)]
    mlist = pd.DataFrame.from_dict(sumlists(mlist['motifs']),
                                   orient='index',
                                   columns=['families'])
    mlist = mlist.sort_values(by='families')
    mlist.plot.bar(legend=False)
    plt.xlabel('Motifs')
    plt.ylabel('Quality-controlled families')
    plt.tight_layout()
    plt.savefig('analysis/qc_fam_motifs.svg')
    plt.close()
                                   
    
def countaa():
    import pandas as pd
    import matplotlib.pyplot as plt
    qual = pd.read_table('data/quality_families.txt')
    qual = qual[qual.total >= 8]
    qual = qual[(qual.experimental_evidence == True) | (qual.perc >= 75)]
    qual = qual['family'].tolist()
    print('Processing absent aas')
    naa = pd.read_table('analysis/absent_aa_norm_family.tsv')
    naa = naa.set_index('family')
    mlist = []
    for row in naa.iterrows():
        n = get_motifs(row)
        if len(n) > 0:
            mlist.append([row[0], n])
            
    mlist = pd.DataFrame(mlist,
                         columns=['family',
                                  'absent_aas'])
    mlist = mlist[mlist.family.isin(qual)]
    mlist = pd.DataFrame.from_dict(sumlists(mlist['absent_aas']),
                                   orient='index',
                                   columns=['families'])
    mlist = mlist.sort_values(by='families')
    mlist.plot.bar(legend=False)
    plt.xlabel('Absent amino acids')
    plt.ylabel('Quality-controlled families')
    plt.tight_layout()
    plt.savefig('analysis/qc_fam_absent_aas.svg')
    plt.close()


def countar():
    import pandas as pd
    import matplotlib.pyplot as plt
    qual = pd.read_table('data/quality_families.txt')
    qual = qual[qual.total >= 8]
    qual = qual[(qual.experimental_evidence == True) | (qual.perc >= 75)]
    qual = qual['family'].tolist()
    print('Processing enriched aas')
    nar = pd.read_table('analysis/rich_aa_norm_family.tsv')
    nar = nar.set_index('family')
    mlist = []
    for row in nar.iterrows():
        n = get_motifs(row)
        if len(n) > 0:
            mlist.append([row[0], n])
            
    mlist = pd.DataFrame(mlist,
                         columns=['family',
                                  'enrich_aas'])
    mlist = mlist[mlist.family.isin(qual)]
    mlist = pd.DataFrame.from_dict(sumlists(mlist['enrich_aas']),
                                   orient='index',
                                   columns=['families'])
    mlist = mlist.sort_values(by='families')
    mlist.plot.bar(legend=False)
    plt.xlabel('Enriched amino acids')
    plt.ylabel('Quality-controlled families')
    plt.tight_layout()
    plt.savefig('analysis/qc_fam_enrich_aas.svg')
    plt.close()


def pa():
    neo, naa, nar = df_generator(merge())
    fixdfs(neo, naa, nar)
    normdfs()
    plot_dfs()
    countfams()
    countaa()
    countar()
        
