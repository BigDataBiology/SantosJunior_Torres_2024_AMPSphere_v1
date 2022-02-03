import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt

# File SPHERE_v.2021-03.levels_assessment.tsv.gz is available in the AMPSphere Zenodo repository
data = pd.read_table('SPHERE_v.2021-03.levels_assessment.tsv.gz')
minsizedfam = [k for k,v in Counter(data['SPHERE_fam level III']).items() if v >= 8]

# File quality_families.txt available in ubuntu@aws.big-data-biology.org:/share/work/Celio/files_for_figures/quality_control
data = pd.read_table('quality_families.txt')
qualfam = data[data.family.isin(minsizedfam)].family.tolist()

# File quality_candidates.txt available in ubuntu@aws.big-data-biology.org:/share/work/Celio/files_for_figures/quality_control
data = pd.read_table('quality_candidates.txt', header=None)
data2 = pd.read_table('high-quality_candidates.txt', header=None)
qualamp = pd.concat([data, data2])[0].tolist()

# Files: families_all.count_core.tsv and amps_all.count_core.tsv
# obtained from analysis_pipe.py run
d = pd.read_table('families_all.count_core.tsv', sep='\t', header='infer')

dfam = pd.DataFrame.from_dict(Counter(d[d.family.isin(minsizedfam)].classification),
                                      orient='index',
                                      columns=['All_families']).T

dqcfam = pd.DataFrame.from_dict(Counter(d[d.family.isin(qualfam)].classification),
                                      orient='index',
                                      columns=['HQ_families']).T


d = pd.read_table('amps_all.count_core.tsv', sep='\t', header='infer')

dfamp = pd.DataFrame.from_dict(Counter(d.classification),
                               orient='index',
                               columns=['All_AMPs']).T
                               
dqcamp = pd.DataFrame.from_dict(Counter(d[d.amp.isin(qualamp)].classification),
                                orient='index',
                                columns=['HQ_AMPs']).T

df = pd.concat([dfamp, dqcamp, dfam, dqcfam]).T

dnew = df * 100 / df.sum()

dnew.T.plot.barh(stacked=True, cmap='Dark2')
plt.xlim(75,100)
plt.xlabel('Percent (%)')
plt.ylabel('Classification')
plt.savefig('core_analysis.svg')
plt.close()

