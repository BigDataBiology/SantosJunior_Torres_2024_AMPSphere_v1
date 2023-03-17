**quality_assessment.tsv**

**Description:**	TSV-file summarizing quality-control results for AMPSphere. 

| **columns** | **description** |
| :---: | :---: |
| AMP | candidate AMP accession code in AMPSphere |
| Antifam | Passed if not matching to AntiFAM database |
| RNAcode | Passed if tested and approved as a protein-coding sequence in RNAcode, however, if not tested, it is referred as such instead of failing |
| metaproteomes | Passed if found in the tested metaproteomes |
| metatranscriptomes | Passed if found in the tested metatranscriptomes |
| Coordinates | Passed if preceeded by a stop codon, inferring that it is not a larger protein fragment |

**MD5 SUM:**	0dc3bc63a5d9c567e5803fecf4fd1ad3

**Size (MBytes):**	1.5968780517578125

**Content sample (first 5 items):**

AMP	Antifam	RNAcode	metaproteomes	metatranscriptomes	Coordinates
AMP10.000_000	Passed	Failed	Passed	Passed	Passed
AMP10.000_001	Passed	Passed	Failed	Passed	Passed
AMP10.000_002	Passed	Failed	Failed	Passed	Passed
AMP10.000_003	Passed	Passed	Passed	Passed	Passed
[...]
