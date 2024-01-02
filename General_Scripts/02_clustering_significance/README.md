# Clusters validation

This analysis consists in the alignment of randomly selected sequences against
the representative sequences for the clusters obtained in the three levels adopted
in the AMPSphere (100%, 100-85%, 100-85-75%).

The set of analysis include:

  (a) selection of random peptide sequences (3000 ptides)
  (b) selection of their representatives, by getting their clusters
  (c) alignment of representative vs. sequence using Smith-Waterman algorithm
  (d) calculate score and e-value using Karlin-Altschman algorithm
  (e) compute significant matches
  (f) plot graphs of distribution and proportion of significant matches by level

To reproduce the analysis:

```bash
python cluster_significance.py
python Fig_clustering_significance.py
jug execute gap-penalty-comparison.py
```

The validation of clusters outputs the following files to the *outputs/* folder:

*clustering_significance.svg*

Scatterplots of Identity (%) versus the Log(E-value)
per alignment per replicate per clustering level.

*output_clustering_significance_levelI.tsv*
*output_clustering_significance_levelII.tsv*
*output_clustering_significance_levelIII.tsv*

These table contains the significance of alignments of randomly selected
sequences against their cluster representatives for the corresponding
clustering levels (100%/85%/75% id). The columns are the following:

|**Column**|**Description**|
| :---: | :---: |
| query	| randomly selected AMP sequence from AMPSphere (not singleton) |
| target | cluster representative |
| identity | amino acids identity without gaps (%) |
| gap_identity | amino acids identity with gaps (%) |
| coverage | coverage of the query sequence |
| score | alignment bit score |
| aln_lenevalue | E-value calculated from the score |
| sig. | significance (* or non-significant/n.s.) |
| family | AMPSphere cluster (Sphere) |


*gap-penalty-comparison.tsv*

Testing different gap penalty parameters

| **Input file** | **Description** |
| :---: | :---: |
| AMPSphere_v.2022-03.faa.gz | AMPSphere peptide sequences |
| SPHERE_v.2022-03.levels_assessment.tsv.gz | Clustering key for the AMPSphere resource |
| AMPsphere_representatives.faa | Representative peptides for each cluster (SPHERE) in the AMPSphere |

