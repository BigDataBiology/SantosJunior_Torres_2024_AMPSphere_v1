**ncbi_missing_lineages.txt**

**Description:**

### The problem

Some lineages are missing in the ProGenomes relase 2. Therefore, we listed
them using the commands:

```
# python3

import pandas as pd

data = pd.read_table('data/AMPSphere_ProGenomes2.tsv.gz', sep='\t', header='infer')
refzero = pd.read_table('data/proGenomes2.1_specI_lineageNCBI.tab', sep='\t', header=None)
refzero = refzero[[0, 7]]
refzero = refzero.rename({0: 'genome', 7: 'source'}, axis=1)
ref = pd.read_table('data/progenomes_samples.tsv', sep='\t', header=None,  names=['specI', 'genome'])
ref = ref.merge(on='genome', right=refzero)
source = []
for x in ref['source'].tolist():
    if x == x: source.append(' '.join(x.split()[1:]))
    else: source.append(np.nan)

ref['source'] = source    
df1 = data.merge(on='genome', right=ref)
df2 = data[~data.GMSC10.isin(df1.GMSC10)]
df2['taxid'] = df2['genome'].apply(lambda x: x.split('.')[0]).astype('int')
taxid_list = df2['taxid'].tolist()
```

By using this taxid_list, we retrieved the lineages
from the [NCBI site](https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi).

Which were later on formated to have just the columns taxid and taxname.

The present file is a TSV-formatted file containing two columns:

 - taxid: NCBI taxonomy identifier 
 - taxname: NCBI taxonomy name taking bionomial nomenclature into account

**MD5 SUM:**	4b45c9d720cabcfc1128f25b43417053

**Size (MBytes):**	0.008998870849609375

**Content sample (first 5 items):**

taxid	taxname
1116126	Escherichia coli 2846750
223967	Methylorubrum populi
1041156	Sinorhizobium meliloti RRl128
1389425	Escherichia coli LAU-EC4
[...]
