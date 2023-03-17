## Code to reproduce the analysis in the AMPSphere manuscript

This collection of Jupyter notebooks bring the code needed to reproduce the research 
in the bioinformatics front from the paper describing the AMPSphere resource.
They are separated by subject and mostly bring a couple of figures and analysis each:

|  **Figure**   |  **Jupyter notebook**  |  **Description**  |
| :---: | :---: | :---: |
|  1A  |  01_Exploring_metadata  |  Map with samples in different colors according the habitat  |
|  1B  |  03_quality_and_homologs_distribution  |  Bar plot and pie chart of AMP quality and mapping per database  |
|  1C  |  04_c_AMP_overlaps  |  Rarefaction curves by environment  |
|  1D  |  02_Graphics_for_habitat_overlap  |  Sankey plots of AMP overlap by habitats  |
|  2A  |  10_gmgc_homologs_analysis  |  Histogram of counts of AMPs versus start of match to full-length proteins as % of the target length  |
|  2B  |  Manual curation  |  Alignment of the gene NAD(P)-dependent dehydrogenase that originates the AMP10.271_016 in different Prevotella genomes  |
|  2CA  |  10_gmgc_homologs_analysis  |  Annotation of AMPs using eggnog mapper – Bar plot showing the functions per COG class  |
|  2CB  |  10_gmgc_homologs_analysis  |  Annotation of AMPs using eggnog mapper – Box plot showing the enrichment in relation to GMGC  |
|  3A  |  11_genome_context  |  Bar plot showing the functions of the more frequent conserved gene neighborhoods involving the AMPs  |
|  3B  |  11_genome_context  |  Bar plot showing the functions involved in the antibiotic resistance in the more frequent conserved gene neighborhoods involving the AMPs  |
|  3C  |  11_genome_context  |  Bar plot showing the functions involving antibiotic synthesis in the more frequent conserved gene neighborhoods involving the AMPs  |
|  3D  |  14_calculate_densities  |  Example for the gene neighborhood of AMP10.015_426  |
|  4A  |  12_clonal_and_accessory_c_AMPs  |  Bar plot showing the proportion of core, shell and accessory AMPs and families with and without high-quality filtering  |
|  4B  |  13_most_represented_taxa  |  Bar plot showing the taxonomic annotation of AMPSphere   |
|  4C  |  16_density_across_taxonomies  |  Box plot showing the AMP density per genus of different phyla  |
|  4D  |  16_density_across_taxonomies  |  Phylogenetic tree of gene found in AMPSphere and showing the AMP density, its associated error and median  |
|  5A  |  01_Exploring_metadata  |  Box plot of AMP density per host vs. non-host associated samples  |
|  5B  |  15_density_across_environments  |  Box plot of AMP density per sample referring to Prevotella copri in different habitats  |
|  5C  |  15_density_across_environments  |  Box plot of average AMP density per species from human oral cavity and guts  |
|  5D  |  15_density_across_environments  |  Box plot of average AMP density per species from soil and plant-associated  |
|  S1A  |  07_c_AMP_features_comparison  |  Distribution of peptide length  |
|  S1B  |  07_c_AMP_features_comparison  |  Distribution of small lateral chain residues (ABCDGNPSTV) in percent  |
|  S1C  |  07_c_AMP_features_comparison  |  Distribution of basic lateral chain residues (HRK) in percent  |
|  S1D  |  07_c_AMP_features_comparison  |  Distribution of pI of peptides  |
|  S1E  |  07_c_AMP_features_comparison  |  Distribution of charge of peptides at pH 7.0  |
|  S1F  |  07_c_AMP_features_comparison  |  Distribution of aliphatic index of peptides  |
|  S1G  |  07_c_AMP_features_comparison  |  Distribution of of instability index of peptides  |
|  S1H  |  07_c_AMP_features_comparison  |  Distribution of Boman index of peptides  |
|  S1I  |  07_c_AMP_features_comparison  |  Distribution of hydrophobic moment of peptides  |
|  S2A  |  03_quality_and_homologs_distribution  |  Bar plot with different peptide quality tests in the y axis and the proportion of AMPSphere in x  |
|  S2B  |  03_quality_and_homologs_distribution  |  Homologs quality by database used in the annotation and the quality test  |
|  S2C  |  04_c_AMP_overlaps  |  Heatmap of the overlap of AMPs across low-level and high-level environments  |
|  S2D  |  05_cAMPs_rarity  |  Line plot of the number of detections versus number of AMPs (genes)  |
|  S3  |  06_c_AMPs_clustering_validation  |  Scatter plots of the identity of the hit against the cluster representative versus the corresponding e-value at the different clustering levels identity cutoffs  |
|  S4  |  01_Exploring_metadata  |  Box plot of AMP density across different high-level habitats  |
|  S5A  |  20_density_across_taxonomies_controlling_quality  |  Box plot of AMP density per genus from different phyla using only quality-controlled AMPs  |
|  S5B  |  21_density_for_host_samples_environments_quality_controlled  |  Box plot of AMP density per sample from different low-level habitats using only quality-controlled AMPs  |
|  S5C  |  21_density_for_host_samples_environments_quality_controlled  |  Box plot of AMP density per sample from host- and non-host-associated environments using only quality-controlled AMPs  |
|  S5D  |  19_density_across_environments_controlling_quality  |  Box plot of AMP density per sample referring to Prevotella copri using only quality-controlled AMPs  |
|  S5E  |  19_density_across_environments_controlling_quality  |  Box plot of average AMP density per species from human oral cavity and guts with only quality-controlled AMPs  |
|  S5F  |  19_density_across_environments_controlling_quality  |  Box plot of average AMP density per species from soil and plant-associated with only quality-controlled AMPs  |


It also brings the Supplementary Tables in the manuscript:

| **Tables**  |  **Jupyter notebook**  |  **Description**  |
| :---: | :---: | :---: |
|  S1  |  01_Exploring_metadata  |  Metadata description associated to the metaproteomes used in this study  |
|  S2  |  09_supplementary_info_about_habitats_and_samples  |  Summarizing statistics about the AMPSphere  |
|  S3  |  10_gmgc_homologs_analysis  |  Orthologs group enrichment in the AMPs that are GMGC homologs  |
|  S4  |  Annotation result  |  KEGG pathway annotation of AMPs with conserved genome contexts  |
|  S5  |  11_genome_context  |  Statistical test of the number of full-length protein families with conserved genome contexts  |
|  S6  |  Annotation result  |  KEGG orthologs to the AMPs with conserved genome contexts  |
|  S7  |  16_density_across_taxonomies  |  AMP density per genus  |
|  S8  |  14_calculate_densities  |  AMP density per species per habitat only in species happening in at least 10 samples per habitat  |
|  S9  |  Manual curation  |  Metadata description associated to the metaproteomes used in this study  |

To open the jupyter notebooks, you will need to type:

```
  $ jupyter notebook
```

