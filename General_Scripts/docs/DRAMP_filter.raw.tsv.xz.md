**DRAMP_filter.raw.tsv**

**Description:**	MMSeqs2 easy-search results using candidate AMPs from AMPSphere against the DRAMP v.3 database.
                        This search was done to filter the singletons matching to a database of validated AMPs, therefore 
                        allowing us to include them in the final AMPSphere resource. The search was raw, meaning that
                        we did filter the results inside the script and the search reported here also contains low-quality
                        hits.

| **columns** | **description** |
| :---: | :---: |
| query | AMP access code from the AMPSphere resource |
| target | AMP access code from the DRAMP v.3 database |
| fident | identity of the alignmet |
| alnlen | alignment length in residues |
| mismatch | number of mismatches in the alignment |
| gapopen | number of gaps openned in the alignment |
| qstart | position of the alignment start in the query sequence |
| qend | position of the alignment stop in the query sequence  |
| tstart | position of the alignment start in the target sequence |
| tend | position of the alignment stop in the target sequence |
| evalue | alignment E-value |
| bits | BitScore calculated from the alignment |

**MD5 SUM:**	51a390b64640eb144244cbeef134273c

**Size (MBytes):**	0.7560043334960938

**Content sample (first 5 items):**

QPGYVPYGAYGGPVPVGCPGGYWARMPVYDPYGNVVGYRGRPRYFCP	DRAMP03794	0.418	38	18	0	1	38	1	32	5.409E-04	31
KTVKVTLLKSISGANRGQRGTVRGLGLKKINHQVELVDTPAVRGMINKVNHLVRVD	DRAMP15043	0.571	56	24	0	1	56	2	57	1.795E-15	64
KTVKVTLLKSISGANRGQRGTVRGLGLKKINHQVELVDTPAVRGMINKVNHLVRVD	DRAMP15039	0.523	56	26	0	1	56	3	58	2.943E-13	58
KTVKVTLLKSISGANRGQRGTVRGLGLKKINHQVELVDTPAVRGMINKVNHLVRVD	DRAMP15046	0.535	54	25	0	3	56	4	57	4.049E-13	58
KTVKVTLLKSISGANRGQRGTVRGLGLKKINHQVELVDTPAVRGMINKVNHLVRVD	DRAMP15042	0.517	56	27	0	1	56	2	57	5.570E-13	57
[...]
