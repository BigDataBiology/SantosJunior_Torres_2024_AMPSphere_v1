# Test of homology

print('Generating protein families dictionary')

seq_dict = dict()
with open(f'{analysis_folder}/AMPsphere_fasta/AMPsphere.faa',
               'rt',
               encoding='utf-8') as db:
    for row in db:
        if row.startswith('>'):
            header = row.strip()
            header = header.replace('>','')
            header = header.split(' | ')
            header = header[0]
        else:
            seq = row.strip()
            seq_dict[header] = seq

pf_dict_I = dict()
pf_dict_II = dict()
pf_dict_III = dict()
pf_dict_IV = dict()
pf_ref_I = dict()
pf_ref_II = dict()
pf_ref_III = dict()
pf_ref_IV = dict()

tab_dic = {'I': 'nr100.sorted.clstr',
           'II': 'nr100-85.sorted.clstr',
           'III': 'nr100-85-75.sorted.clstr',
           'IV': 'nr100-85-75-50.sorted.clstr'}

pfs_rel = {'I': [pfs_I, pf_dict_I, pf_ref_I],
           'II': [pfs_II, pf_dict_II, pf_ref_II],
           'III': [pfs_III, pf_dict_III, pf_ref_III],
           'IV': [pfs_IV, pf_dict_IV, pf_ref_IV]}

for s in ['I', 'II', 'III', 'IV']:
    data = pd.read_table(f'{analysis_folder}/clustering/{tab_dic[s]}',
                         sep='\t',
                         header=None)
    for i in range(len(data[0])):
        printProgressBar(i, len(data[0]), 'sequences from sphere level-I')
        if data.iloc[i][4] == '*':
            pfs_rel[s][2][data.iloc[i][0]] = data.iloc[i][1]
        else:
            pfs_rel[s][1][data.iloc[i][1]] = data.iloc[i][0]
    pfs_rel[s][0] = sorted(pfs_rel[s][1].values())
    pfs_rel[s][0] = [key for key, group in groupby(pfs_rel[s][0]) if len(list(group)) > 5]
    print('Global pair-wise alignment\n')
    print(f'\t\tPreparing fasta and references\n')
    total = len(pfs_rel[s][0])
    i = 1
    nfasta1 = []
    for protfam in pfs_rel[s][0]:
        printProgressBar(i, total, f' -- processed {i} protein families')
        seqlist = [key for key, value in pfs_rel[s][1].items() if value == protfam]
        n = [(key, seq_dict[key]) for key in seqlist]
        nfasta1 += n
        i += 1
    nfasta1 = set(nfasta1)
    print(f'\t\tReplicate 1 of SPHERE_{s} with 1000 randomly sampled sequences\n')
    fasta_sub = random.sample(nfasta1, 1000)  # subsampling fasta for replicate 1
    aln_test(fasta_sub, '1', s)
    print(f'\t\tReplicate 2 of SPHERE_{s}  with 1000 randomly sampled sequences\n')
    fasta_sub = random.sample(nfasta1, 1000)  # subsampling fasta for replicate 1
    aln_test(fasta_sub, '2', s)
    print(f'\t\tReplicate 3 of SPHERE_{s} with 1000 randomly sampled sequences\n')
    fasta_sub = random.sample(nfasta1, 1000)  # subsampling fasta for replicate 1
    aln_test(fasta_sub, '3', s)

