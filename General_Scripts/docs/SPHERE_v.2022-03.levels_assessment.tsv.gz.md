**SPHERE_v.2022-03.levels_assessment.tsv**

**Description:**	Clustering key of candidate AMPs from AMPSphere. The TSV-file contains
                        per AMP the different clusters they belong using an hierarchical system
                        where level I represents a clustering at 100% of identity, level II 
                        corresponds to the clustering of those representatives at 85% of identity,
                        and the third level represents a cluster of the remainant representatives
                        at 75%. To that, it was used a reduced alphabet to improve the clustering
                        power of distant homologs as described by [Murphy et al. (2000)](https://doi.org/10.1093/protein/13.3.149).

| **columns** | **description** |
| :---: | :---: |
| AMP accession | AMP unique accession from AMPSphere |
| evaluation vs. representative | representative sequence ('*') or the alignment with the representative (qstart:qstop:tstart:tstop/identity) |
| SPHERE_fam level I | cluster identifier at level I |
| SPHERE_fam level II | cluster identifier at level II |
| SPHERE_fam level III | cluster identifier at level III |

**MD5 SUM:**	a1c45d76db598fa861f8725380f87c74

**Size (MBytes):**	12.306156158447266

**Content sample (first 5 items):**

AMP accession	evaluation vs. representative	SPHERE_fam level I	SPHERE_fam level II	SPHERE_fam level III
AMP10.000_000	*	SPHERE-I.000_766	SPHERE-II.000_794	SPHERE-III.001_493
AMP10.000_001	*	SPHERE-I.002_463	SPHERE-II.000_195	SPHERE-III.000_082
AMP10.000_002	*	SPHERE-I.248_008	SPHERE-II.008_664	SPHERE-III.009_856
AMP10.000_003	1:34:1:36/94.12%	SPHERE-I.000_444	SPHERE-II.000_080	SPHERE-III.000_200
[...]
