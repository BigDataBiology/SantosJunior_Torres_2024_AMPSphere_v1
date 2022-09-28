**quality_families.txt**

**Description:**	TSV-file summarizing the quality control results for AMPs at the family level.
                        To this, we consider the clusters obtained at the third level (100-85-75% of identity).

| **columns** | **description** |
| :---: | :---: |
| family | cluster accession in AMPSphere |
| experimental_evidence | presence in metaproteomes or metatranscriptomes tested |
| quality_candidates | number of AMPs in the cluster passing all quality tests¹ |
| total | number of AMPs composing the referred cluster |
| perc | percent of candidates passing all quality tests¹ |

*NOTE:* ¹ - quality tests here refers to the AntiFam, RNAcode and coordinates test.
In case the experimental evidence is obtained the entire cluster is considered high-quality, 
analogously if the 'perc' >= 75%.

**MD5 SUM:**	a75f9072cb61abf735d234318ba347b8

**Size (MBytes):**	0.49497222900390625

**Content sample (first 5 items):**

family	experimental_evidence	quality_candidates	total	perc
SPHERE-III.000_000	True	543	1790	30.34
SPHERE-III.000_001	True	842	872	96.56
SPHERE-III.000_002	True	527	546	96.52
SPHERE-III.000_003	True	405	497	81.49
[...]
