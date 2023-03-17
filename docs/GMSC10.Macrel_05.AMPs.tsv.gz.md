**GMSC10.Macrel_05.AMPs.tsv**

**Description:**	Candidate antimicrobial peptides predicted with Macrel from the 
                        predicted small ORFs from GMSC. It is a TSV-file containing
                        only positive sequences.

| **columns** | **description** |
| :---: | :---: |
| Access | gene identifier in the GMSC resource |
| Sequence | predicted peptide sequence | 
| AMP_family | classification of peptide according to macrel¹ |
| AMP_probability | probability of being a true positive using macrel's model |
| Hemolytic | class for hemolysis activity (Hemo/NonHemo) |
| Hemolytic_probability | probability of being hemolytic according to macrel's model |

¹ AMP family classified into anionic/cationic and cysteine-containing/linear coded as:
  A/C and D/L, respectively, followed by a P (for peptide), _e.g._, a cationic
  cysteine-containing peptide would be CDP

**MD5 SUM:**	b3c3bcf345731624c070644902e92110

**Size (MBytes):**	263.29407787323

**Content sample (first 5 items):**

Access	Sequence	AMP_family	AMP_probability	Hemolytic	Hemolytic_probability
GMSC10.SMORF.000_000_000_338	RGGVSKNSHMDGSGPDQNPGHKAGMMFYQNFSQVMRSYMST	CLP	0.644	NonHemo	0.198
GMSC10.SMORF.000_000_000_949	KLRKILKSMFNNYCKTFKDVPPGNMFR	CDP	0.614	Hemo	0.772
GMSC10.SMORF.000_000_001_095	MKGYKHVNKKPKLIPQIYYAVKPLQ	CLP	0.574	Hemo	0.9109999999999999
GMSC10.SMORF.000_000_002_097	CFTSKNCYKTGNFSCVLKQKDPDFLKRFSAIWFF	CDP	0.6729999999999999	NonHemo	0.436
[...]
