## Relation of inputs and outputs per analysis 

Here is the relation of input and output files used and produced by each core analysis 
in the AMPSphere workflow. The outputs here mentioned are those actively used in the 
analysis making to the paper and available in the
[GitHub repository](https://github.com/BigDataBiology/AMPSphere_manuscript/tree/main/Manuscript_Analysis)
as jupyter notebooks.

The description of files can be found in the links available.

**01_resource_generation**

*Input*
        
 - [alph.json](alph.json.xz.md)
 - [AMPsphere_representatives.faa](AMPsphere_representatives.faa.xz.md)
 - [AMPSphere_v.2022-03.faa.gz](AMPSphere_v.2022-03.faa.gz.md)
 - [cmap.json](cmap.json.xz.md)
 - [combinedVsearch3.50p.20n.norm.opt.0.042.filteredBlackList.map](combinedVsearch3.50p.20n.norm.opt.0.042.filteredBlackList.map.xz.md)
 - [DRAMP.fa](DRAMP.fa.xz.md)
 - [dramp_anno.py](dramp_anno.py.xz.md)
 - [DRAMP_filter.raw.tsv](DRAMP_filter.raw.tsv.xz.md)
 - [GMSC10.Macrel_05.AMPs.tsv.gz](GMSC10.Macrel_05.AMPs.tsv.gz.md)
 - [GMSC10.ProGenomes2.coords.txt.gz](GMSC10.ProGenomes2.coords.txt.gz.md)
 - [GMSC10.proGenomesv2_05.AMPs.rename.tsv.gz](GMSC10.proGenomesv2_05.AMPs.rename.tsv.gz.md)
 - [gmsc_genes.fna.xz](gmsc_genes.fna.xz.md)
 - [metadata.tsv](metadata.tsv.xz.md)
 - [proGenomes2.1_specI_clustering.tab](proGenomes2.1_specI_clustering.tab.xz.md)
 - [proGenomes2.1_specI_lineageNCBI.tab](proGenomes2.1_specI_lineageNCBI.tab.xz.md)
 - [SPHERE_v.2022-03.levels_assessment.tsv.gz](SPHERE_v.2022-03.levels_assessment.tsv.gz.md)
	
*Output*
	
 - [SPHERE_v.2022-03.levels_assessment.tsv.gz](SPHERE_v.2022-03.levels_assessment.tsv.gz.md)
 - [AMPSphere_v.2022-03.fna.xz](AMPSphere_v.2022-03.fna.xz.md)
 - [AMPsphere_representatives.faa](AMPsphere_representatives.faa.xz.md)
 - [AMPSphere_v.2022-03.faa.gz](AMPSphere_v.2022-03.faa.gz.md)
 - [ampsphere_v2022-03.features.tsv.gz](ampsphere_v2022-03.features.tsv.gz.md)
 - [dramp_v3.features.tsv.gz](dramp_v3.features.tsv.gz.md)
 - [macrel_trainpos.features.tsv.gz](macrel_trainpos.features.tsv.gz.md)

**02_clustering_significance**

*Input*
        
 - [AMPsphere_representatives.faa](AMPsphere_representatives.faa.xz.md)
 - [AMPSphere_v.2022-03.faa.gz](AMPSphere_v.2022-03.faa.gz.md)
 - [SPHERE_v.2022-03.levels_assessment.tsv.gz](SPHERE_v.2022-03.levels_assessment.tsv.gz.md)
	
*Output*
	
 - [output_clustering_significance_levelI.tsv.gz](output_clustering_significance_levelI.tsv.gz.md)
 - [output_clustering_significance_levelII.tsv.gz](output_clustering_significance_levelII.tsv.gz.md)
 - [output_clustering_significance_levelIII.tsv.gz](output_clustering_significance_levelIII.tsv.gz.md)

**03_metadata_analysis**

*Input*
        
 - [AMPSphere_v.2022-03.fna.xz](AMPSphere_v.2022-03.fna.xz.md)
 - [AMPSphere_v.2022-03.origin_samples.tsv.gz](AMPSphere_v.2022-03.origin_samples.tsv.gz.md)
 - [ar122_metadata_r202.tsv](ar122_metadata_r202.tsv.xz.md)
 - [bac120_metadata_r202.tsv](bac120_metadata_r202.tsv.xz.md)
 - [general_envo_names.tsv](general_envo_names.tsv.xz.md)
 - [GMSC10.metag_smorfs.rename.txt.xz](GMSC10.metag_smorfs.rename.txt.xz.md)
 - [GMSC10.ProGenomes2.coords.txt.gz](GMSC10.ProGenomes2.coords.txt.gz.md)
 - [metadata.tsv](metadata.tsv.xz.md)
 - [mmseqs2.lca_taxonomy.full.tsv.xz](mmseqs2.lca_taxonomy.full.tsv.xz.md)
 - [ncbi_missing_lineages.txt](ncbi_missing_lineages.txt.xz.md)
 - [proGenomes2.1_specI_clustering.tab](proGenomes2.1_specI_clustering.tab.xz.md)
 - [proGenomes2.1_specI_lineageNCBI.tab](proGenomes2.1_specI_lineageNCBI.tab.xz.md)	
        
*Output*
	
 - [complete_amps_associated_taxonomy.tsv.gz](complete_amps_associated_taxonomy.tsv.gz.md)
 - [gmsc_amp_genes_envohr_source.tsv.gz](gmsc_amp_genes_envohr_source.tsv.gz.md)
 - [metadata.tsv](metadata.tsv.xz.md)
 - [reduced_metadata.tsv](reduced_metadata.tsv.xz.md)

**04_rnacode_amp_families**

*Input*
        
 - [AMPSphere_v.2022-03.fna.xz](AMPSphere_v.2022-03.fna.xz.md)
 - [SPHERE_v.2022-03.levels_assessment.tsv.gz](SPHERE_v.2022-03.levels_assessment.tsv.gz.md)

*Output*

 - [RNAcode_out_wlfam.tsv](RNAcode_out_wlfam.tsv.xz.md)

**05_quality**

*Input*
        
 - [AMPSphere_v.2022-03.faa.gz](AMPSphere_v.2022-03.faa.gz.md)
 - [AMPSphere_v.2022-03.fna.xz](AMPSphere_v.2022-03.fna.xz.md)
 - [antifam_100.tsv.gz](antifam_100.tsv.gz.md)
 - [coordinates_test_passed.tsv.gz](coordinates_test_passed.tsv.gz.md)
 - [DRAMP_annotation.raw.tsv](DRAMP_annotation.raw.tsv.xz.md)
 - [gmsc_amp_genes_envohr_source.tsv.gz](gmsc_amp_genes_envohr_source.tsv.gz.md)
 - [metaproteome_100.tsv.gz](metaproteome_100.tsv.gz.md)
 - [metaT_100.tsv.gz](metaT_100.tsv.gz.md)
 - [RNAcode_out_wlfam.tsv](RNAcode_out_wlfam.tsv.xz.md)
 - [SPHERE_v.2022-03.levels_assessment.tsv.gz](SPHERE_v.2022-03.levels_assessment.tsv.gz.md)

*Output*

 - [quality_assessment.tsv](quality_assessment.tsv.xz.md).gz
 - [high_[quality_candidates.txt](quality_candidates.txt.xz.md)](high_[quality_candidates.txt](quality_candidates.txt.xz.md).xz.md).gz
 - [quality_candidates.txt](quality_candidates.txt.xz.md).gz
 - [quality_families.txt](quality_families.txt.xz.md).gz

**06_homologs_and_overlap**

*Input*
        
 - [AMPSphere_v.2022-03.faa.gz](AMPSphere_v.2022-03.faa.gz.md)
 - [databases_homology](databases_homology.tar.xz.md)
 - [general_envo_names.tsv](general_envo_names.tsv.xz.md)
 - [gmsc_amp_genes_envohr_source.tsv.gz](gmsc_amp_genes_envohr_source.tsv.gz.md)
 - [high_[quality_candidates.txt](quality_candidates.txt.xz.md)](high_[quality_candidates.txt](quality_candidates.txt.xz.md).xz.md)
 - [metadata.tsv](metadata.tsv.xz.md)
 - [quality_assessment.tsv](quality_assessment.tsv.xz.md)
 - [quality_candidates.txt](quality_candidates.txt.xz.md)

*Output*

 - [adjust_significant_function.csv.xz](adjust_significant_function.csv.xz.md)
 - [amp_COG.tsv.xz](amp_COG.tsv.xz.md)
 - [dramp_candidates.txt.gz](dramp_candidates.txt.gz.md)
 - [gmgc_candidates.txt.gz](gmgc_candidates.txt.gz.md)
 - [SmProt_candidates.txt.gz](SmProt_candidates.txt.gz.md)
 - [starPepDB_candidates.txt.gz](starPepDB_candidates.txt.gz.md)
 - [STsORFs_candidates.txt.gz](STsORFs_candidates.txt.gz.md)
 - [result_gmgc.m8.xz](result_gmgc.m8.xz.md)

**07_taxonomy_core_analysis**

*Input*
        
 - [37185.fasta](37185.fasta.xz.md)
 - [bac120_r202.tre](bac120_r202.tre.xz.md)
 - [bac120_taxonomy_r202.tsv](bac120_taxonomy_r202.tsv.xz.md)
 - [bps-per-sample-per-taxon.tsv.xz](bps-per-sample-per-taxon.tsv.xz.md)
 - [bps-per-taxon.tsv.xz](bps-per-taxon.tsv.xz.md)
 - [complete_amps_associated_taxonomy.tsv.gz](complete_amps_associated_taxonomy.tsv.gz.md)
 - [gmsc_amp_genes_envohr_source.tsv.gz](gmsc_amp_genes_envohr_source.tsv.gz.md)
 - [quality_candidates.txt](quality_candidates.txt.xz.md)
 - [high_quality_candidates.txt](high_quality_candidates.txt.xz.md)
 - [pgenomes_samples.tsv](pgenomes_samples.tsv.xz.md)
 - [prevotella_species_list.tsv](prevotella_species_list.tsv.xz.md)
 - [quality_candidates.txt](quality_candidates.txt.xz.md)
 - [quality_families.txt](quality_families.txt.xz.md)
 - [SPHERE_v.2022-03.levels_assessment.tsv.gz](SPHERE_v.2022-03.levels_assessment.tsv.gz.md)

*Output*

 - [amps_all.count_core.tsv.gz](amps_all.count_core.tsv.gz.md)
 - [families_all.count_core.tsv.gz](families_all.count_core.tsv.gz.md)

**08_amps_in_progenomes_ANI_core**

*Input*

 - [complete_gmsc_pgenomes_metag.tsv](complete_gmsc_pgenomes_metag.tsv.xz.md)
 - [proGenomes2.1_specI_clustering.tab](proGenomes2.1_specI_clustering.tab.xz.md)
 - [proGenomes2.1_specI_lineageNCBI.tab](proGenomes2.1_specI_lineageNCBI.tab.xz.md)
 - [SPHERE_v.2022-03.levels_assessment.tsv.gz](SPHERE_v.2022-03.levels_assessment.tsv.gz.md)

*Output*
	
 - [amp_results_ani.tsv](amp_results_ani.tsv.xz.md)
 - [summary_output_core_prots.tsv](summary_output_core_prots.tsv.xz.md)

**09_AMPs_from_Mpneumoniae**

*Input*

 - [fam_red.tsv](fam_red.tsv.xz.md)
 - [genomes_list](genomes_list.xz.md)
 - [genotyping_Mpneumoniae_Diaz2017_journal.pone.0174701.tsv](genotyping_Mpneumoniae_Diaz2017_journal.pone.0174701.tsv.xz.md)
 - [mgpa_ref.fa](mgpa_ref.fa.xz.md)
 - [progenomes_cluster1969.tsv](progenomes_cluster1969.tsv.xz.md)

*Output*
	
 - [sample_genes](sample_genes.tar.xz.md)
 - [tree_alignment](tree_alignment.tar.xz.md)
 - [type1.txt](type1.txt.xz.md)
 - [type2.txt](type2.txt.xz.md)
 - [type3.txt](type3.txt.xz.md)

**10_rarity**

*Input*

 - [uniques.tar.gz](uniques.tar.gz)
 - [envo.pkl](envo.pkl.md)
 - [headers.tsv.xz](headers.tsv.xz.md)

*Output*
	
 - [unique_map_cAMPs_high_habs.tsv.xz](unique_map_cAMPs_high_habs.tsv.xz.md)
 - [unique_map_cAMPs_general_habs.tsv.xz](unique_map_cAMPs_general_habs.tsv.xz.md)
 - [unique_map_cAMPs_nsamples.tsv.xz](unique_map_cAMPs_nsamples.tsv.xz.md)

**11_host_non_host_amps**

*Input*

 - [general_envo_names.tsv](general_envo_names.tsv.xz.md)
 - [gmsc_amp_genes_envohr_source.tsv.gz](gmsc_amp_genes_envohr_source.tsv.gz.md)
 - [metadata.tsv](metadata.tsv.xz.md)
 - [samples-min500k-assembly-prodigal-stats.tsv](samples-min500k-assembly-prodigal-stats.tsv.xz.md)

*Output*
	
 - [host_factor.svg](host_factor.svg.md)
 - [amp-density-per-habitat.svg](amp-density-per-habitat.svg.md)
 - [mannwhitneyu_test_mammalguts.tsv](mannwhitneyu_test_mammalguts.tsv.md)

**12_genome_context**

*Output*
	
 - [AMPs_hit_genomes.txt.xz](AMPs_hit_genomes.txt.xz.md)
 - [Genomic_context_frequency_KEGG_pathways.tsv.xz](Genomic_context_frequency_KEGG_pathways.tsv.xz.md)
 - [hits.folded.tab.xz](hits.folded.tab.xz.md)
 - [neigh_VCS_higher_0.9_3_or_more_hits.tab.xz](neigh_VCS_higher_0.9_3_or_more_hits.tab.xz.md)
 - [permutation_tests.tab.xz](permutation_tests.tab.xz.md)
 - [strong_hits.txt.xz](strong_hits.txt.xz.md)


**13_transmissibility**

*Input*

 - [complete_amps_associated_taxonomy.tsv.gz](complete_amps_associated_taxonomy.tsv.gz.md)
 - [metadata.tsv.xz](metadata.tsv.xz.md)

*Output*
	
 - [species_amp_richness_crossenvironment.tsv.gz](species_amp_richness_crossenvironment.tsv.gz.md)

**14_BGI_peptides**

*Input*
        
 - [AMPSphere_v.2022-03.faa.gz](AMPSphere_v.2022-03.faa.gz.md)
 - [BGI_peptides.faa](BGI_peptides.faa.xz.md)

*Output*

 - [patents_blast.tsv](patents_blast.tsv.xz.md)
 - [patents_diamond.tsv](patents_diamond.tsv.xz.md)
 - [patents_mmseqs.tsv](patents_mmseqs.tsv.xz.md)
 - [matches/](matches.tar.gz.md)
    
**15_select_candidates_for_synthesis**

*Input*

 - [AMPs_per_genomes.tsv](AMPs_per_genomes.tsv.xz.md)
 - [ampsphere_v2022-03.features.tsv.gz](ampsphere_v2022-03.features.tsv.gz.md)
 - [AMPSphere_v.2022-03.solsyn_rules.tsv.gz](AMPSphere_v.2022-03.solsyn_rules.tsv.gz.md)
 - [ensemble_predictions.tsv](ensemble_predictions.tsv.xz.md)
 - [gmsc_amp_genes_envohr_source.tsv.gz](gmsc_amp_genes_envohr_source.tsv.gz.md)
 - [quality_assessment.tsv](quality_assessment.tsv.xz.md)
 - [result_dramp.m8](result_dramp.m8.xz.md)
 - [result_SmProt.m8](result_SmProt.m8.xz.md)
 - [result_starPepDB.m8](result_starPepDB.m8.xz.md)
 - [SPHERE_v.2022-03.levels_assessment.tsv.gz](SPHERE_v.2022-03.levels_assessment.tsv.gz.md)

*Output*
	
 - [AMPSphere_features.tsv](AMPSphere_features.tsv.md)
 - [AMPSphere_fmt_involved.tsv](AMPSphere_fmt_involved.tsv.md)
 - [AMPSphere_homologs.tsv](AMPSphere_homologs.tsv.md)
 - [AMPSphere_quality_control.tsv](AMPSphere_quality_control.tsv.md)
 - [AMPSphere_selected_candidates_coprediction.tsv](AMPSphere_selected_candidates_coprediction.tsv.md)
 - [AMPSphere_sequences.tsv](AMPSphere_sequences.tsv.md)
 - [AMPSphere_synthesis_solubility.tsv](AMPSphere_synthesis_solubility.tsv.md)
 - [AMPSphere_taxonomy_environment.tsv](AMPSphere_taxonomy_environment.tsv.md)
 - [selected_candidates.tsv](selected_candidates.tsv.md)
 - [selected_peptides.fasta](selected_peptides.fasta.md)

