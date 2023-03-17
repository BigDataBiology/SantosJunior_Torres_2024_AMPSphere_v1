**patents_diamond.tsv**

**Description:**	Comparison of patented sequences from [Ma et al. (2021)](https://www.nature.com/articles/s41587-022-01226-0),
                        available in the patent linked [here](https://app.dimensions.ai/details/patent/WO-2022037681-A1),
                        against AMPSphere peptides by using [Diamond software](https://github.com/bbuchfink/diamond). The output
                        consists of a TSV-file with the following non-named columns:

| **column position** | **description** |
| 1 | query sequence from the patent |
| 2 | query amino acids sequence |
| 3 | query sequence length in residues |
| 4 | target sequence from AMPSphere |
| 5 | target amino acids sequence |
| 6 | target sequence length |
| 7 | identity 0-100% |
| 8 | alignment length |
| 9 | E-value |
| 10 | query coverage |
| 11 | target coverage |

**MD5 SUM:**	00f9d382b4728655c641de6bad78f4c9

**Size (MBytes):**	0.010532379150390625

**Content sample (first 5 items):**

AMP2	GKGIKFVGEEIRRKSGKSAGAK	22	AMP10.132_201	VAAKIRSLRKPEPYKGKGIKFVGEQLRRKAGKSAGAK	37	86.4	22	4.25e-07	100	59.5
AMP2	GKGIKFVGEEIRRKSGKSAGAK	22	AMP10.238_699	LGQVAAKIRSFRKPEPYKGKGIKFVGEEIRRKAGKSAGKK	40	90.9	22	4.49e-07	100	55.0
AMP2	GKGIKFVGEEIRRKSGKSAGAK	22	AMP10.624_462	ICAKIRSFRKPEPYKGKGVKFVGEIIRRKSGKSAGAK	37	90.9	22	6.04e-07	100	59.5
AMP2	GKGIKFVGEEIRRKSGKSAGAK	22	AMP10.340_042	VCAKIRSFRKPEPYKGKGILFVGEQIRRKSGKSAGAK	37	90.9	22	8.58e-07	100	59.5
AMP2	GKGIKFVGEEIRRKSGKSAGAK	22	AMP10.634_509	VCAKIRSFRKPEPYKGKGILFLGEQIRRKSGKSAGAK	37	86.4	22	2.46e-06	100	59.5
[...]
