**GMSC10.proGenomesv2_05.AMPs.rename.tsv**

**Description:**	AMP prediction using Macrel in the ProGenomes v2 samples.
                        TSV-file containing the raw names of genes and their properties
                        during prediction.

| **columns** | **description** |
| :---: | :---: |
| Access | small ORF identifier, containing taxid, bioproject, contig, and arbitrary ORF identifier |
| Sequence | predicted amino acids sequence |
| AMP_family | classification of peptide according to macrel¹ |
| AMP_probability | probability of being a true positive using macrel's model |
| Hemolytic | class for hemolysis activity (Hemo/NonHemo) |
| Hemolytic_probability | probability of being hemolytic according to macrel's model |

¹ AMP family classified into anionic/cationic and cysteine-containing/linear coded as:
  A/C and D/L, respectively, followed by a P (for peptide), _e.g._, a cationic
  cysteine-containing peptide would be CDP

**MD5 SUM:**	e66b19b78fd75acb45f217ea73b20463

**Size (MBytes):**	0.8681764602661133

**Content sample (first 5 items):**

# Prediction from macrel v0.5.0
Access	Sequence	AMP_family	AMP_probability	Hemolytic	Hemolytic_probability
1262903.PRJEB687.FR903869_4_28	RGGVSKNSHMDGSGPDQNPGHKAGMMFYQNFSQVMRSYMST	CLP	0.644	NonHemo	0.198
344610.PRJNA15639.AAKB02000001_4120	KLRKILKSMFNNYCKTFKDVPPGNMFR	CDP	0.614	Hemo	0.772
439842.PRJNA19457.ABAK02000001_173	MKGYKHVNKKPKLIPQIYYAVKPLQ	CLP	0.574	Hemo	0.911
[...]
