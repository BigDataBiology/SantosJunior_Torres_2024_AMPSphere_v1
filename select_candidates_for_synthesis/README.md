**Solubility Rules:**

Solubility rules are based on recommendations of [SolyPep](http://bioserv.rpbs.univ-paris-diderot.fr/services/SolyPep/)

```
	(1) A peptide cannot start or end with a charged residue
	(2) A peptide cannot contain more than 1 P or G in the sequence
	(3) A peptide should not be composed more than 25% by a single residue
	(4) A peptide should not have more than 45% of charged/hydrophobic residues
	(5) A peptide should not have more than 75% of gel-prone residues
	(6) A peptide should not have net charge above 1
```

**Synthesis Rules:**

Synthesis rules are based on empirical observations from the [PepFun](https://github.com/rochoa85/PepFun) group

```
	(1) A peptide should not have prohibited motifs: 2 consecutive prolines, DG or DP, or N/Q at amino-terminal position
	(2) Charged residues should not constitute more than 50% of the residue
	(3) No residues that are oxidized (CMW)
```


**Inputs:**

  - AMPSphere_v.2022-03.faa.gz:	Fasta of AMPSphere peptides
  
  - SPHERE_v.2022-03.levels_assessment.tsv.gz:	Table containing the association of AMPs and their families in AMPSphere

  - high_quality_candidates.txt:	AMP candidates with experimental evidence of translation/transcription and passing in all 
  other tests (Antifam, Terminal placement, RNACode)

  - quality_candidates.txt:	AMP candidates without experimental evidence of translation/transcription and passing in all 
  other tests (Antifam, Terminal placement, RNACode)


**Outputs:**

Results are organized as three tables, they both contain the same columns:

| **Column** | **Description** |
| :---: | :---: |
| accession | AMPSphere peptide accession code |
| sol_rule_1 | Solubility rule 1 |
| sol_rule_2 | Solubility rule 2 |
| sol_rule_3 | Solubility rule 3 |
| sol_rule_4 | Solubility rule 4 |
| sol_rule_5 | Solubility rule 5 |
| sol_rule_6 | Solubility rule 6 |
| syn_rule_1 | Synthesis rule 1 |
| syn_rule_2 | Synthesis rule 2 |
| syn_rule_3 | Synthesis rule 3 |
| fam | SPHERE family at level III |

Tables:

  - AMPSphere_v.2022-03.solsyn_rules.tsv.gz:	contains all peptides from AMPSphere and their assessment regarding solubility and synthesis
  
  - AMPSphere_v.2022-03.selected_candidates.tsv:	contains only quality controlled peptides from AMPSphere and their assessment regarding sol. and synthesis
  
  - selected_candidates.tsv:	table contains the candidates suggested for synthesis given their superior results, this table contains an additional column
  FScore which consists in the geometric mean of the percent of True results for synthesis and solubility, if it was above 60, then the candidate was selected.

