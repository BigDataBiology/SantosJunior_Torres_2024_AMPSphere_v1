# Motifs.py

Internal script
Run:

```
python3 uscripts/motifs.py
```

### : Description :

It generates a series of tables describing the motifs distribution
by families and high-quality families. Also generates tables bringing
the distribution of these motifs by families. 

### : Inputs :

1. **data/SPHERE_v.2022-03.levels_assessment.tsv.gz:**

Table from AMPSphere v2022-03 version. It consists in a table with the clustering 
information in three different levels (100%, 100-85%, 100-85-75%) after conversion 
using an 8-letter alphabet. Its columns are:

    AMP accession - accession in AMPSphere v2022-03 resource
    evaluation vs. representative - summary of clustering procedure
    SPHERE_fam level I - cluster accession at 1st level (100% identity)
    SPHERE_fam level II - cluster accession at 2nd level (100-85% identity)
    SPHERE_fam level III - cluster accession at 3rd level (100-85-75% identity)

2. **data/quality_families.txt.gz:**

Table generated during the quality assessment. It contains the number of clusters
passing all the tests and which percent it represents for the cluster, also
contains the information about experimental evidence of at least 1 candidate in the
cluster. Its columns are:

    family - AMP cluster accession at 3rd level (100-85-75% identity)
    experimental_evidence - True or False('*')
    quality_candidates - number of peptides passing in all quality tests ('**')
    total - total number of peptides in the cluster
    perc - percent of candidates passing in all tests

'*'  regarding to the identification of the AMP in metaproteomes or metatranscriptomes
   at least 1 candidate needed to become True

'**' simulatenously in RNAcode, AntiFam, and Terminal placement

3. **data/AMPSphere_v.2022-03.annotation.tsv.gz:**

Table generated during annotation of motifs. The general regex expressions were searched
using the function __str__.contains() and results are shown by presence/absence.
Its columns are:

    id - accession in AMPSphere v2022-03 resource
    AA_rich - most frequent residue (abundance %)
    AA_absent - absent amino acids in the peptide
    sol.1 - does not start or end with a charged residue
    sol.2 - does not contain more than 1 P or G in the sequence
    sol.3 - is not composed more than 25% by a single residue
    sol.4 - no more than 45% of charged/hydrophobic residues
    sol.5 - no more than 75% of gel-prone residues
    sol.6 - net charge <= 1
    syn.1 - no forbiden motifs: 2 consecutive prolines,
            DG or DP, or N/Q at amino-terminal position
    syn.2 - charged residues <= 50% of peptide
    syn.3 - no residues that are oxidized (CMW)
    motif_match - str, contains the list of motifs per AMP ('*')
    
'*' motifs searched were available in [Ruhanen et al. (2014)](https://pubmed.ncbi.nlm.nih.gov/24478765/)
    and [Huan et al. (2021)](https://pubmed.ncbi.nlm.nih.gov/33178164/). To a complete list
    see the following `data/db_motif.tsv.gz`.
    
### : Outputs :

1. **motif_seqs_per_family_counts.tsv**

Table bringing the number of AMPs per motif (in columns) per family (in rows).
Only clusters with at least 8 peptides were kept for this analysis.

2. **motif_seqs_per_family_norm.tsv**

Table bringing the percent of the family contains a given motif (in columns)
per family (in rows). Only clusters with at least 8 peptides were kept for
this analysis.

3. **motifs_listed_per_family.tsv.gz**

Table containing two columns `family` and `motifs`. The families consist
of clusters of at least 8 peptides. The motifs presented in this table
happened in at least 75% of the peptides from a given cluster.

4. **counts_motifs_per_family.tsv.gz**

Table with two columns `Motifs` and `High-quality families`. The column
`High-quality families` brings the number of clusters of at least 8 peptides
passing in all quality tests which have at least 75% of their members with
that motif.

5. **qc_fam_motifs.svg**

Horizontal bar chart made with the table `counts_motifs_per_family.tsv.gz`, in 
which the `x axis` is the number of high-quality families and `y axis` is the
motif.


