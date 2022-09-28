## Search for AMP motifs

This analysis used the database compiled from the studies
of Ruhuanen et al. and Huan et al. and uses as inputs the files:

 - [AMPSphere_v.2022-03.faa.gz](../data_folder/docs/AMPSphere_v.2022-03.faa.gz.md)
 - [db_motif.tsv](../data_folder/docs/db_motif.tsv.xz.md)
 - [environment_classification.py](../data_folder/docs/environment_classification.py.xz.md)
 - [gmsc_amp_genes_envohr_source.tsv.gz](../data_folder/docs/gmsc_amp_genes_envohr_source.tsv.gz.md)
 - [quality_families.txt](../data_folder/docs/quality_families.txt.xz.md)
 - [SPHERE_v.2022-03.levels_assessment.tsv.gz](../data_folder/docs/SPHERE_v.2022-03.levels_assessment.tsv.gz.md)

Using ReGex rules and the python 3.9, we obtained the following outputs:

### Tables

Brings amino acids enrichment per sequence as proportions for
AMPs on their clusters - `rich_aa_seqs_per_family.tsv`
and later normalized by the number of sequences presenting those
enrichments - `rich_aa_norm_family.tsv`

Brings the motifs present in the sequences - `motif_seqs_per_family.tsv`
and families - `motifs_norm_family.tsv`

Brings the amino acids grouped by their features per sequence
and environment - `groupped_amino_acids_by_seq_by_env.tsv`

Brings a complete table with AMPs and their motifs and other
important features involving amino acids distribution and
depletion - `AMPSphere_v.2022-03.annotation.tsv`

Brings the distribution of amino acids in the AMPs 
found in different environments - `amino_acids_distribution_by_seq_by_env.tsv`

Also brings the depleted amino acids in the sequences
per family - `absent_aa_seqs_per_family.tsv` and 
normalizes this by the size of the families - `absent_aa_norm_family.tsv`

### Figures

amino_acids_groupedI_composition.svg
motif_top_sps.svg
qc_fam_absent_aas.svg
qc_fam_enrich_aas.svg
qc_fam_motifs.svg
motif_envo.svg
motif_mammal_guts.svg
motif_seqs.svg
enrich_aa_seqs.svg
amino_acids_composition.svg
amino_acids_groupedI_composition.svg
absent_aa_seqs.svg

### Folders:

**clustermaps/:**
