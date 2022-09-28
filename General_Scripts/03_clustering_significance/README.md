# Clusters validation

This analysis consists in the alignment of randomly selected sequences against
the representative sequences for the clusters obtained in the three levels adopted
in the AMPSphere (100%, 100-85%, 100-85-75%). This is important to validate the 
reducted alphabet we used to increase the coalescence of clusters.

The set of analysis include:

  (a) selection of random peptide sequences (1000 in 3 replicates)
  (b) selection of their representatives, by getting their clusters
  (c) alignment of representative vs. sequence using Smith-Waterman algorithm
  (d) calculate score and e-value using Karlin-Altschman algorithm
  (e) compute significant matches
  (f) plot graphs of distribution and proportion of significant matches by level

The validation was done using two different approaches, one with the complete 
amino acids alphabet and another with the reduced alphabet (redalph).

To reproduce the analysis:

```
# to reproduce the validation of clusters:
    $ python3 main.py
```

The validation of clusters outputs the following files to the *analysis/* folder:

*clustering_significance.svg*

Scatterplots of Identity (%) versus the Log(E-value)
per alignment per replicate per clustering level.

*coverage_distribution.svg*

Density plot of coverage showing the distribution of
coverages by alignments by clustering levels.

*identity_distribution.svg*

Distribution of Identity (%) of alignments per clustering
level.

*redalph_identity_distribution.svg*

Distribution of Identity (%) of alignments per clustering
level by using the reduced alphabet.

*redalph_significance_across_clustering_levels.svg*

Bar plot showing the proportion of significant alignments
per clustering level by using the reduced alphabet.

*significance_across_clustering_levels.svg*

Bar plot showing the proportion of significant alignments
per clustering level.

*output_clustering_significance_levelI.tsv*

This table contains the significance of alignments of randomly 
selected sequences against their cluster representatives for the
clustering level I (100% id). The columns are the following:
    
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
| replicate | round of replication of random sampling |

*output_clustering_significance_levelII.tsv*

This table contains the significance of alignments of randomly 
selected sequences against their cluster representatives for the
clustering level II (100-85% id). The columns are the following:
    
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
| replicate | round of replication of random sampling |

*output_clustering_significance_levelIII.tsv*

This table contains the significance of alignments of randomly 
selected sequences against their cluster representatives for the
clustering level III (100-85-75% id). The columns are the following:
    
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
| replicate | round of replication of random sampling |

*output_clustering_redalph_significance_levelI.tsv*

This table contains the significance of alignments of randomly 
selected sequences against their cluster representatives after
a reduction of the amino acids alphabet for the
clustering level I (100% id). The columns are the following:
    
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
| replicate | round of replication of random sampling |

*output_clustering_redalph_significance_levelII.tsv*

This table contains the significance of alignments of randomly 
selected sequences against their cluster representatives after
a reduction of the amino acids alphabet for the
clustering level II (100-85% id). The columns are the following:
    
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
| replicate | round of replication of random sampling |

*output_clustering_redalph_significance_levelIII.tsv*

This table contains the significance of alignments of randomly 
selected sequences against their cluster representatives after
a reduction of the amino acids alphabet for the
clustering level III (100-85-75% id). The columns are the following:
    
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
| replicate | round of replication of random sampling |
    
The scripts in *utils/* folder are explained below:

	- clustering_analysis.py:**Function.** Perform clustering validation of randomly selected sequences
	                         **Inputs.** SPHERE_v.2022-03.levels_assessment.tsv.gz, AMPsphere_representatives.faa, AMPSphere_v.2022-03.faa.gz
	                         **Outputs.** output_clustering_significance_levelI.tsv, output_clustering_significance_levelII.tsv, output_clustering_significance_levelIII.tsv

	- clustering_redalph.py: **Function.** Perform clustering validation of randomly selected sequences after amino acids alphabet reduction
	                         **Inputs.** SPHERE_v.2022-03.levels_assessment.tsv.gz, AMPsphere_representatives.faa, AMPSphere_v.2022-03.faa.gz
	                         **Outputs.** output_clustering_redalph_significance_levelI.tsv, output_clustering_redalph_significance_levelII.tsv, output_clustering_redalph_significance_levelIII.tsv

	- displotter.py:         **Function.** Plot the distribution of features in the clustering validation procedures
	                         **Inputs.** output_clustering_significance_levelI.tsv, output_clustering_significance_levelII.tsv, output_clustering_significance_levelIII.tsv
	                         **Outputs.** clustering_significance.svg

	- plot_output.py:        **Function.** Plot proportion distribution of clustering significance at different levels
	                         **Inputs.** output_clustering_significance_levelI.tsv, output_clustering_significance_levelII.tsv, output_clustering_significance_levelIII.tsv
	                         **Outputs.** significance_across_clustering_levels.svg, identity_distribution.svg, coverage_distribution.svg

	- plot_redalph.py:       **Function.** Plot proportion distribution of clustering significance using reduced amino acids alphabet at different levels
	                         **Inputs.** output_clustering_redalph_significance_levelI.tsv, output_clustering_redalph_significance_levelII.tsv, output_clustering_redalph_significance_levelIII.tsv
	                         **Outputs.** redalph_significance_across_clustering_levels.svg, redalph_identity_distribution.svg, redalph_coverage_distribution.svg

Detailed inputs list:

| **Input file** | **Description** |
| :---: | :---: |
| AMPSphere_v.2022-03.faa.gz | AMPSphere peptide sequences |
| SPHERE_v.2022-03.levels_assessment.tsv.gz | Clustering key for the AMPSphere resource |
| AMPsphere_representatives.faa | Representative peptides for each cluster (SPHERE) in the AMPSphere |

