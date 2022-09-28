    
from itertools import chain
amps_found = set(list(chain.from_iterable(amp_data.sequences.dropna().str.split(',').tolist())))
from Bio import SeqIO
selected_seqs = []
for i in selected_seqs:
    if i in amps_found: print(i)

amp_data = amp_data[['genome', 'species',
                     'sample', 'fmt.id',
                     'ORFs', 'smORFs',
                     '#_amps']]

merge_donor = merge_donor.rename({'sample.pre': 'sample'}, axis=1)
merge_donor = merge_donor[['genome', 'species', '#_amps', 'outcome']]
merge_donor = amp_data.merge(on=['species', 'sample', 'fmt.id'], right=merge_donor)
merge_donor =merge_donor[['genome', 'species', '#_amps', 'outcome']]
merge_donor = merge_donor.sort_values(by='species')
checkdf(merge_donor)
'''
coexistence	vs.	species retained: U=65942.5, p=1.5130708867745961e-05
engraftment conspecific	vs.	species lost: U=3118.0, p=0.048744059438622256
engraftment conspecific	vs.	species retained: U=10289.5, p=0.0004922256392564286
rejection conspecific	vs.	species retained: U=13989.5, p=0.03304375341763645
species lost	vs.	engraftment conspecific: U=2162.0, p=0.048744059438622256
species retained	vs.	coexistence: U=44473.5, p=1.5130708867745961e-05
species retained	vs.	engraftment conspecific: U=6190.5, p=0.0004922256392564286
species retained	vs.	rejection conspecific: U=10730.5, p=0.03304375341763645
'''


merge_recipient = result.drop('sample', axis=1)
merge_recipient = merge_recipient.rename({'sample.post': 'sample'}, axis=1)
merge_recipient = merge_recipient[['genome', 'species', '#_amps', 'outcome']]
merge_recipient = amp_data.merge(on=['species', 'sample', 'fmt.id', 'timepoint.fmt'], right=merge_recipient)
merge_recipient = merge_recipient[['genome', 'species', '#_amps', 'outcome']]
merge_recipient = merge_recipient.sort_values(by='species')
checkdf(merge_recipient)
'''
coexistence	vs.	influx novel: U=11653.5, p=0.01992394479997522
coexistence	vs.	species retained: U=42779.0, p=0.02322739148169393
engraftment novel	vs.	influx novel: U=5473.5, p=0.04424207516480239
influx novel	vs.	coexistence: U=7696.5, p=0.01992394479997522
influx novel	vs.	engraftment novel: U=3771.5, p=0.04424207516480239
species retained	vs.	coexistence: U=34171.0, p=0.02322739148169393
'''

recipient = []
for record in merge_recipient.groupby(['species', 'outcome']):
    recipient.append([record[0][0], record[0][1], record[1]['#_amps'].tolist()])

recipient = pd.DataFrame(recipient, columns=['species', 'outcome', 'amps'])
recipient['type'] = 'recipient'

donor = []
for record in merge_donor.groupby(['species']):
    donor.append([record[0], '-', record[1]['#_amps'].tolist()])

donor = pd.DataFrame(donor, columns=['species', 'outcome', 'amps'])
donor['type'] = 'donor'

df = pd.concat([recipient, donor])
df = df.sort_values(by=['species', 'outcome'])
df.to_csv('fmt_amps_per_genome_and_outcome.tsv', sep='\t', header=True, index=None)

df = df[df.amps.apply(lambda x: len(x) >= 10)]
df = df[df.species.isin([k for k,v in Counter(df.species).items() if v > 1])]
for record in df.groupby('species'):
    for j in record[1].outcome:
        if j != '-':
            u, p = mu(record[1][record[1].outcome == '-']['amps'].tolist()[0], record[1][record[1].outcome == j]['amps'].tolist()[0])
            if p < 0.05:
                print(record[0])
                print(f'donor vs. recipient -- {j}: u={u}, p={p:.4f}')
'''
ref_mOTU_v2_0947
donor vs. recipient -- engraftment novel: u=75.0, p=0.0175
ref_mOTU_v2_2358
donor vs. recipient -- coexistence: u=522.0, p=0.0232
ref_mOTU_v2_4266
donor vs. recipient -- engraftment conspecific: u=321.5, p=0.0335
ref_mOTU_v2_4266
donor vs. recipient -- engraftment novel: u=280.0, p=0.0189
'''

sp_class = amp_data.merge(on=['species', 'sample', 'fmt.id'], right=result)

sp_class = sp_class[['genome', 'species', 'sample', 'fmt.id', 'ORFs', 'smORFs', '#_amps',
                     'fmt.type', 'timepoint.fmt', 'sample.pre', 'sample.post', 'outcome']]

sp_class.to_csv('results_per_outcome.tsv', sep='\t', header=True, index=None)

checkdf(sp_class)
'''
coexistence	vs.	engraftment novel: U=46930.5, p=2.850903027610371e-05
coexistence	vs.	rejection novel: U=30580.0, p=0.009998956977171531
engraftment conspecific	vs.	engraftment novel: U=7075.5, p=0.00922102899348932
engraftment novel	vs.	coexistence: U=31080.5, p=2.850903027610371e-05
engraftment novel	vs.	engraftment conspecific: U=4689.5, p=0.00922102899348932
engraftment novel	vs.	rejection conspecific: U=4202.0, p=0.00013377122516182823
engraftment novel	vs.	species lost: U=2067.0, p=7.735530980036972e-05
rejection conspecific	vs.	engraftment novel: U=7744.0, p=0.00013377122516182823
rejection conspecific	vs.	rejection novel: U=5063.5, p=0.0039165354382420655
rejection novel	vs.	coexistence: U=22864.0, p=0.009998956977171531
rejection novel	vs.	rejection conspecific: U=3120.5, p=0.0039165354382420655
rejection novel	vs.	species lost: U=1558.5, p=0.0013226191562018127
species lost	vs.	engraftment novel: U=4630.0, p=7.735530980036972e-05
species lost	vs.	rejection novel: U=3029.5, p=0.0013226191562018127
'''

x = t['#_amps'].apply(lambda x: Counter(x))
newdf = pd.DataFrame()
for idx in x.index:
    n = [(idx, k, v*100/sum(x.loc[idx].values())) for k, v in x.loc[idx].items()]
    n = pd.DataFrame(n, columns=['outcome', 'AMPs', 'pct'])
    newdf = pd.concat([newdf, n])

newdf = newdf.pivot_table(index='outcome', columns='AMPs', values='pct').fillna(0)
newdf.to_csv('percent_distribution_per_outcome.tsv', sep='\t', header=True, index=True)

x = sp_class[sp_class.species.isin([k for k, v in Counter(sp_class.species).items() if v >= 10])]

x[['species',
   'fmt.id',
   'ORFs',
   'smORFs',
   '#_amps',
   'timepoint.fmt',
   'outcome']].to_csv('outcome.tsv',
                      sep='\t',
                      header=True,
                      index=None)

for record in x[['species', 'outcome', '#_amps']].groupby('species'):
    print(record[0])
    t = record[1][['#_amps', 'outcome']]
    t = t.groupby('outcome').agg(list)
    t = t.loc[t['#_amps'].apply(lambda x: len(x)) > 3]
    test_groups(t)

'''
meta_mOTU_v2_5364
engraftment novel	vs.	coexistence: U=72.0, p=0.007279949258093863

meta_mOTU_v2_6452
engraftment conspecific	vs.	rejection novel: U=25.0, p=0.028254027550050152
engraftment novel	vs.	coexistence: U=17.5, p=0.03688842570704987
engraftment novel	vs.	rejection novel: U=30.0, p=0.022950538964290012
rejection conspecific	vs.	rejection novel: U=25.0, p=0.028254027550050152

meta_mOTU_v2_6557
coexistence	vs.	rejection conspecific: U=156.0, p=0.014987840869210893

meta_mOTU_v2_7093
coexistence	vs.	species lost: U=58.0, p=0.02177316855360182

ref_mOTU_v2_1301
rejection novel	vs.	engraftment novel: U=15.5, p=0.03247332032606962

ref_mOTU_v2_1416
rejection conspecific	vs.	coexistence: U=376.0, p=0.0011880898717912594

'''

ax = sns.boxplot(x='outcome', y='#_amps', data=x)
ax = sns.swarmplot(x="outcome", y="#_amps", data=x, color=".25")
ax.set_xticklabels(ax.get_xticklabels(),rotation = 30)
plt.tight_layout()
plt.savefig('outcome_analysis.svg')

test = test[test.outcome.isin([k for k, v in Counter(test.outcome).items() if v >= 10])]
ax = sns.boxplot(x="outcome", y="#_amps", data=test)
ax = sns.swarmplot(x="outcome", y="#_amps", data=test, color=".25", s=0.5)
ax.set_xticklabels(ax.get_xticklabels(),rotation = 30)
plt.tight_layout()
plt.savefig('outcome_anlaysis_all.svg')

