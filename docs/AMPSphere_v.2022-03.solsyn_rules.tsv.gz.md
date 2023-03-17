**AMPSphere_v.2022-03.solsyn_rules.tsv**

**Description:**	gives information about the criteria the AMPs
                        passed or failed for solubility and synthesis.

Solubility rules are based on recommendations of [SolyPep](http://bioserv.rpbs.univ-paris-diderot.fr/services/SolyPep/).

Synthesis rules are based on empirical observations from the [PepFun](https://github.com/rochoa85/PepFun) group.

| **columns** | **description** |
| :---: | :---: |
| accession | access code in the AMPSphere resource |
| sol_rule_1 | True if not start or end with a charged residue |
| sol_rule_2 | True if not contain more than 1 P or G in the sequence |
| sol_rule_3 | True if not be composed more than 25% by a single residue |
| sol_rule_4 | True if not have more than 45% of charged/hydrophobic residues |
| sol_rule_5 | True if not have more than 75% of gel-prone residues |
| sol_rule_6 | True if not have net charge above 1 |
| syn_rule_1 | True if not have prohibited motifs: 2 consecutive prolines, DG or DP, or N/Q at amino-terminal position |
| syn_rule_2 | True if charged residues should not constitute more than 50% of the residue |
| syn_rule_3 | True if no residues that are oxidized (CMW) |
| fam | candidate AMP cluster at 3rd level (100-85-75% identity) |

**MD5 SUM:**	76273ed8b07af3725692001e0e44fcdc

**Size (MBytes):**	5.954895973205566

**Content sample (first 5 items):**

accession	sol_rule_1	sol_rule_2	sol_rule_3	sol_rule_4	sol_rule_5	sol_rule_6	syn_rule_1	syn_rule_2	syn_rule_3	fam
AMP10.000_000	False	False	False	False	True	False	True	False	False	SPHERE-III.001_493
AMP10.000_001	True	False	True	False	True	False	True	False	False	SPHERE-III.000_082
AMP10.000_002	False	False	True	False	True	False	True	False	False	SPHERE-III.009_856
AMP10.000_003	False	False	True	False	True	False	True	False	False	SPHERE-III.000_200
[...]
