**bps-per-sample-per-taxon.tsv**

**Description:**	TSV-file containing the number of base pairs assembled per taxonomic groups per sample
                        for the GMSC resource from the metagenome datasets. Taxonomy was obtained using 
                        the approach of MMSeqs2 with the GTDB database.

| **column** | **description** |
| :---: | :---: |
| sample | metagenome sample identifier |
| taxid | Taxonomy ID as in NCBI database |
| level | Taxonomy rank |
| name | Taxonomy group name (binomial when available) |
| nbps | Number of base pairs assembled for a given taxonomy group in a given sample |

**MD5 SUM:**	1f5595877d53ae01f57b5fd5155288bd

**Size (MBytes):**	1202.1895141601562

**Content sample (first 5 items):**

sample	taxid	level	name	nbps
Karasov_2018_arabidopsis_NextMet1	0	no rank	unclassified	159873
Karasov_2018_arabidopsis_NextMet1	1	no rank	root	40460
Karasov_2018_arabidopsis_NextMet1	2	superkingdom	Bacteria	426285
Karasov_2018_arabidopsis_NextMet1	3	phylum	Proteobacteria	34237
[...]
