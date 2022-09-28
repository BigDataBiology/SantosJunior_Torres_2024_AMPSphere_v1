# AMPSphere Generation

This analysis integrates all steps needed to generate AMPSphere resources.
It intends to generate a set of tables that besides useful for analysis,
also are used in the data management.

The set of analysis include:

	(a) selection of non-singleton peptides
	(b) recovery of singletons matching DRAMP
	(c) reduction of alphabet
	(d) hierarchical clustering
	(e) calculation of features
	(f) calculate families alignment
	(g) calculate families tree and ascII representation
	(h) calculate HMM profile and logo
	(i) calculate helical wheels
	(j) add progenomes info and metadata info
	(k) generation of files to deposit in Zenodo

To reproduce the analysis:

```
# to reproduce AMPSphere:
    $ python3 main.py
```

Three fastas are outputted to the *analysis/* folder:

| **Fasta file** | **Description** |
| :---: | :---: |
| AMPsphere_representatives.faa | This fasta contains the peptide sequences for the representatives for each cluster |
| AMPSphere_v.2022-03.faa.gz | This fasta contains the peptide sequences for all AMPs forming AMPSphere |
| AMPSphere_v.2022-03.fna.xz | This fasta contains the nucleotide sequences for the genes found coding for AMPs in AMPSphere |

The folders outputted in *analysis/* folder are:

	clustering - this folder contains files generated during the hierarchical clustering with CDHIT

	families - this folder contains other 6 folders built with info from clusters at level III with at least 8 peptides:

		fastas - contains fasta files of peptide sequences from each cluster
		
		aln - contains the alignment built with Mafft for each cluster
		
		hmm - contains the Hidden Markov Models (HMM) built from the alignments
		
		hmm_logo - contains HMM logos built for each model
		
		tree_nwk - contains the Newick files for trees built using FastTree2 for each cluster

		tree_fig - contains the ASCII images for each tree
		
		seqlogo - contains sequence logos for each alignment
		
	helical_wheels - it contains the image files in .svg format for the helical wheels built for each AMP in AMPSphere

AMPSphere generation procedure outputs the following tables to the *analysis/* folder:

*AMPsphere_GMSC_correspondence.tsv.gz*

This table contains GMSC genes and the AMPs ending up the AMPSphere.
The columns are the following:
    
|**Column**|**Description**|
| :---: | :---: |
| accession | Antimicrobial peptide accession in AMPSphere |
| sequence | peptide sequence for the AMP in AMPSphere |
| length | length of the AMP peptide in residues |
| n_of_genes | number of genes coding for a given AMP in AMPSphere |
| genes | genes from GMSC10 resource coding for AMPs in AMPSphere |

*AMPsphere_metaG_annotation.tsv.gz*

This table contains GMSC genes and the AMPs ending up the AMPSphere.
The columns are the following:
    
|**Column**|**Description**|
| :---: | :---: |
| accession | Antimicrobial peptide accession in AMPSphere |
| gmsc | GMSC10 gene accession |
| sample | metagenome sample |
| contig | assembled contig id |
| start | gene start position in the contig |
| stop | gene stop position in the contig |
| strand | orientation of the gene, frame (1/-1) |
| fid | ordinary identity of the gene |
| partial | partial gene (10/01) or complete (00) |
| start_type | start codon |
| rbs_motif | ribosome binding site motif |
| rbs_spacer | ribosome binding site spacer |
| gc_cont | Guanine/Cytosine percent content in the gene  |
    
*AMPsphere_proGenomes_correspondence.tsv.gz*

This table contains AMPs from ProGenomes with their genome of origin and specI cluster.
The columns are the following:
    
|**Column**|**Description**|
| :---: | :---: |
| accession | Antimicrobial peptide accession in AMPSphere |
| Sequence | peptide sequence |
| genome | genome accession from which the AMP is originating from |
| specI | species cluster |

*AMPsphere_species_pGenomes.tsv.gz*

This table contains AMPs from ProGenomes with their genome of origin and specI cluster,
and NCBI lineage details.
The columns are the following:
    
|**Column**|**Description**|
| :---: | :---: |
| accession | Antimicrobial peptide accession in AMPSphere |
| Sequence | peptide sequence |
| genome | genome accession from which the AMP is originating from |
| specI | species cluster |
| kingdom | superkingdom for the genome belongs to |
| phylum | phylum for the genome belongs to |
| class | class for the genome belongs to |
| order | order for the genome belongs to |
| family | family for the genome belongs to |
| genus | genus for the genome belongs to |
| species | species for the genome belongs to |

*AMPSphere_v.2022-03.features.tsv.gz*

This table contains AMPs from ProGenomes with their genome of origin and specI cluster,
and NCBI lineage details.
The columns are the following:
    
|**Column**|**Description**|
| :---: | :---: |
| accession | Antimicrobial peptide accession in AMPSphere |
| tinyAA | percent of amino acids with tiny lateral chain (A + C + G + S + T) |
| smallAA | percent of amino acids with small lateral chain (A + B + C + D + G + N + P + S + T + V) |
| aliphaticAA | percent of amino acids with aliphatic lateral chain (A + I + L + V) |
| aromaticAA | percent of amino acids with aromatic lateral chain (F + H + W + Y) |
| nonpolarAA | percent of amino acids with non-polar lateral chain (A + C + F + G + I + L + M + P + V + W + Y) |
| polarAA | percent of amino acids with polar lateral chain (D + E + H + K + N + Q + R + S + T + Z) |
| chargedAA | percent of amino acids with charged lateral chain (B + D + E + H + K + R + Z) |
| basicAA | percent of amino acids with basic lateral chain (H + K + R) |
| acidicAA | percent of amino acids with acidic lateral chain (B + D + E + Z) |
| charge | net charge of the peptide at pH 7.0 |
| pI | isoelectric point calculated from the lateral chain pKas |
| aindex | aliphatic index calculated as relative volume occupied by aliphatic side chains (A, V, I, and L) |
| instaindex | instability index which predicts stability of a protein based on its amino acid composition, if lower than 40 is predicted as stable |
| boman | Boman index: numeric vector with the potential peptide-interaction index of the amino acid sequence |
| hydrophobicity | hydrophobicity measured with Einsenberg scale |
| hmoment | Hidrophobic moment measuring of the amphiphilicity perpendicular to the axis of a potential alpha helix |
| SA.Group1.residue0 | CTD-D at first residue of the surface accessibility amino acids group 1 |
| SA.Group2.residue0 | CTD-D at first residue of the surface accessibility amino acids group 2 |
| SA.Group3.residue0 | CTD-D at first residue of the surface accessibility amino acids group 3 |
| HB.Group1.residue0 | CTD-D at first residue of the transfer energy amino acids group 1 |
| HB.Group2.residue0 | CTD-D at first residue of the transfer energy amino acids group 2 |
| HB.Group3.residue0 | CTD-D at first residue of the transfer energy amino acids group 3 |

*AMPSphere_v.2022-03.hosts.tsv.gz*

This table contains the number of genes per AMP found per host.
The columns are the following:
    
|**Column**|**Description**|
| :---: | :---: |
| accession | Antimicrobial peptide accession in AMPSphere |
| host_common_name |  popular name of the host species |
| host_scientific_name | scientific name of the host species |
| host_tax_id | NCBI taxonomy identification |
| counts | number of genes counted for the AMP in a given host |

*AMPSphere_v.2022-03.locations.tsv.gz*

This table contains the number of genes per AMP found in each geographic location.
The columns are the following:
    
|**Column**|**Description**|
| :---: | :---: |
| accession | Antimicrobial peptide accession in AMPSphere |
| geographic_location |  geographic location, such as countries (China, e.g.) and oceans (Atlantic, e.g.) |
| counts | number of genes counted for the AMP in a given location |

*AMPSphere_v.2022-03.microontology.tsv.gz*

This table contains the number of genes per AMP found in each microontology term.
The columns are the following:
    
|**Column**|**Description**|
| :---: | :---: |
| accession | Antimicrobial peptide accession in AMPSphere |
| microontology |  microontology term to classify a given sample |
| counts | number of genes counted for the AMP in a given microontology term |

*AMPSphere_v.2022-03.origin_samples.tsv.gz*

This table contains the genes and samples from metagenomes and ProGenomes for each AMP in AMPSphere.
The columns are the following:
    
|**Column**|**Description**|
| :---: | :---: |
| accession | Antimicrobial peptide accession in AMPSphere |
| genes | list of GMSC10 genes codifying for the AMP |
| progenomes | genomes containing a given AMP from ProGenomes |
| metagenomes | metagenome samples in which a given AMP was found |

*AMPSphere_v.2022-03.species.tsv.gz*

This table contains AMPs from ProGenomes with their genome of origin and specI cluster.
The columns are the following:
    
|**Column**|**Description**|
| :---: | :---: |
| accession | Antimicrobial peptide accession in AMPSphere |
| genome | genome accession from which the AMP is originating from |
| specI | species cluster |

*DRAMP_anno_AMPSphere_v.2022-03.parsed.tsv.gz*

Blast-like table contaning singleton with homology to DRAMP
The columns are the following:
    
|**Column**|**Description**|
| :---: | :---: |
| query | Antimicrobial peptide accession in AMPSphere |
| target | DRAMP homolog |
| fident | percent identity |
| alnlen | alignment length |
| mismatch | number of mismatches |
| gapopen | number of gaps opened in the alignemnt |
| qstart | Position in the query sequence for the start of the alignment |
| qend | Position in the query sequence for the end of the alignment |
| tstart | Position in the target sequence for the start of the alignment |
| tend | Position in the target sequence for the end of the alignment |
| evalue | E-value for the alignment, significance |
| bits | Bit score |
| description | Name of the peptide and origin |
| activity | Description of antimicrobial activity |
| phase_of_testing | Phase of the clinical test |
| reference | Works published describing the peptide |

*DRAMP_annotation.raw.tsv*

Blast-like table contaning singleton with homology to DRAMP, raw as it 
was generated by MMSeqs
The columns are the following:
    
|**Column**|**Description**|
| :---: | :---: |
| query | Antimicrobial peptide accession in AMPSphere |
| target | DRAMP homolog |
| fident | percent identity |
| alnlen | alignment length |
| mismatch | number of mismatches |
| gapopen | number of gaps opened in the alignemnt |
| qstart | Position in the query sequence for the start of the alignment |
| qend | Position in the query sequence for the end of the alignment |
| tstart | Position in the target sequence for the start of the alignment |
| tend | Position in the target sequence for the end of the alignment |
| evalue | E-value for the alignment, significance |
| bits | Bit score |

*DRAMP_filter.parsed.tsv*

Blast-like table contaning singleton with homology to DRAMP after the parsing
procedure where hits were filtered by significance and identity.
The columns are the following:
    
|**Column**|**Description**|
| :---: | :---: |
| query | Antimicrobial peptide accession in AMPSphere |
| target | DRAMP homolog |
| fident | percent identity |
| alnlen | alignment length |
| mismatch | number of mismatches |
| gapopen | number of gaps opened in the alignemnt |
| qstart | Position in the query sequence for the start of the alignment |
| qend | Position in the query sequence for the end of the alignment |
| tstart | Position in the target sequence for the start of the alignment |
| tend | Position in the target sequence for the end of the alignment |
| evalue | E-value for the alignment, significance |
| bits | Bit score |

*SPHERE_v.2022-03.levels_assessment.tsv.gz*

Table contaning peptide clusters association at different levels
The columns are the following:
    
|**Column**|**Description**|
| :---: | :---: |
| AMP accession | Antimicrobial peptide accession in AMPSphere |
| evaluation vs. representative | Clustering information |
| SPHERE_fam level I | Cluster at 100% of identity |
| SPHERE_fam level II | Cluster at 100-85% of identity |
| SPHERE_fam level III | Cluster at 100-85-75% of identity |
    
The scripts in *utils/* folder are explained below:

	- timeout_input.py:      **Function.** Contains the function of timeout when requesting an user input
	                         **Inputs.** None
	                         **Outputs.** None

	- utils.py:              **Function.** Contains the functions for call different softwares (cdhit, crev, mmseqs, hmmbuild)
	                         **Inputs.** None
	                         **Outputs.** None

	- phylogen.py:           **Function.** Contains functions for the tree computation using FastTree2 and its representation as ASCII
	                         **Inputs.** None
	                         **Outputs.** None

	- dramp_anno.py:         **Function.** DRAMP complementar resource. Dictionary structured resource.
	                         **Inputs.** None
	                         **Outputs.** None

	- database_features.py:  **Function.** Contains information and constants used in the AMP features calculation.
	                         **Inputs.** None
	                         **Outputs.** None

	- hmm.py:                **Function.** Brings functions for HMM profile building and logo rendering
	                         **Inputs.** None
	                         **Outputs.** None

	- singletons_handle.py:  **Function.** Selects non-singleton AMPs and recover the singletons matching to DRAMP
	                         **Inputs.** DRAMP_GCP2019.fasta, GMSC10.Macrel_05.AMPs.tsv.gz
	                         **Outputs.** AMPsphere_representatives.faa, AMPSphere_v.2022-03.faa.gz, clustering/,
	                         AMPsphere_GMSC_correspondence.tsv.gz, SPHERE_v.2022-03.levels_assessment.tsv.gz, DRAMP_filter.parsed.tsv

	- features.py:           **Function.** Computes peptide features
	                         **Inputs.** AMPs from AMPSphere (AMPSphere_v.2022-03.faa.gz)
	                         **Outputs.** AMPSphere_v.2022-03.features.tsv.gz

	- annotation.py:         **Function.** Annotates AMPSphere using DRAMP database by a MMSeqs2 search
	                         **Inputs.** Database of AMPs (DRAMP_GCP2019.fasta), AMPs from AMPSphere (AMPSphere_v.2022-03.faa.gz)
	                         **Outputs.** Results of the search (DRAMP_annotation.raw.tsv) and parsed results (DRAMP_anno_AMPSphere_v.2022-03.parsed.tsv.gz)

	- progenomes_amps.py:    **Function.** Associate progenomes info to the AMP and state their origins
	                         **Inputs.** GMSC10.proGenomesv2_05.AMPs.rename.tsv.gz, AMPsphere_GMSC_correspondence.tsv.gz, pgenomes_samples.tsv,
	                         combinedVsearch3.50p.20n.norm.opt.0.042.filteredBlackList.map, AMPsphere_proGenomes_correspondence.tsv.gz, proGenomes2.1_specI_lineageNCBI.tab
	                         **Outputs.** AMPsphere_proGenomes_correspondence.tsv.gz, AMPsphere_species_pGenomes.tsv.gz

	- cluster_analysis.py:   **Function.** Analize families of AMPs generating HMM profiles, logo, trees, tree representations, alignments, fasta file per family
	                         **Inputs.** AMPSphere_v.2022-03.faa.gz
	                         **Outputs.** families/fastas, families/aln, families/hmm, families/hmm_logo, families/tree_nwk, families/tree_fig

	- helical.py:            **Function.** Compute helical wheels for each AMP in AMPSphere, generate a .svg per sequence
	                         **Inputs.** AMPSphere_v.2022-03.faa.gz
	                         **Outputs.** helical_wheels/

	- metaG.py:              **Function.** Associate metadata info to the metagenome AMPs
	                         **Inputs.** AMPsphere_GMSC_correspondence.tsv.gz, GMSC10.metag_smorfs.rename.txt.xz, AMPsphere_proGenomes_correspondence.tsv.gz, metadata.tsv
	                         **Outputs.** AMPsphere_metaG_annotation.tsv.gz, AMPSphere_v.2022-03.hosts.tsv.gz, AMPSphere_v.2022-03.locations.tsv.gz,
	                         AMPSphere_v.2022-03.microontology.tsv.gz, AMPSphere_v.2022-03.origin_samples.tsv.gz, AMPSphere_v.2022-03.species.tsv.gz

	- gmsc_genes.py:         **Function.** Generate genes fasta for AMPSphere
	                         **Inputs.** AMPsphere_GMSC_correspondence.tsv.gz, gmsc_genes.fna.xz
	                         **Outputs.** AMPSphere_v.2022-03.fna.xz

	- seqlogo.py:            **Function.** Generate sequence logos for each peptide alignment in families/aln folder
	                         **Inputs.** Alignment files in families/aln/ automatically retrievied and listed
	                         **Outputs.** PDF files in families/seqlogo/ folder, containing the sequence logo for each corresponding alignment.

Detailed inputs list:

| **Input file** | **Description** |
| :---: | :---: |
| alph.json | JSON file containing color shades for the HMMLogo package |
| cmap.json | JSON file containing color shades for the HMMLogo package |
| combinedVsearch3.50p.20n.norm.opt.0.042.filteredBlackList.map | Map of high-quality genomes made by Askaberg from EMBL for further filtering of ProGenomes resource |
| dramp_anno.py | Manually curated resource of DRAMP with the annotations for each AMP in the database as a dictionary |
| DRAMP.fa | [DRAMP v.3](DRAMP.cpu-bioinfor.org/) |
| DRAMP_filter.raw.tsv | Pre-computed resource for AMPSphere using the old MMSeqs version, it is a blast-like table |
| GMSC10.Macrel_05.AMPs.tsv.gz | JSON file containing color shades for the HMMLogo package |
| GMSC10.ProGenomes2.coords.txt.gz | Internal GMSC resource. Table with coordinates for all predicted smORFs in the genomes from ProGenomes2 |
| GMSC10.proGenomesv2_05.AMPs.rename.tsv.gz | Internal GMSC resource. Table with the AMPs from ProGenomes database, it contains the results from Macrel v.0.5: 'Access\|Sequence\|AMP_family\|AMP_probability\|Hemolytic\|Hemolytic_probability' |
| GMSC10.metag_smorfs.rename.txt.xz | Internal GMSC resource. Table with all coordinates and predictions of GMSC genes from metagenomes, with accession codes to GMSC database |
| gmsc_genes.fna.xz | GMSC10 filtered genes for AMPSphere, pre-computed resource made available by Luis Pedro Coelho |
| metadata.tsv | Table containing the metadata for the metagenome samples used to predict smORFs, currently it is allocated in the [public repository](https://github.com/BigDataBiology/global_data/tree/master/freeze.v2) |
| proGenomes2.1_specI_clustering.tab | ProGenomes v.2 specI and samples table, available in the resource through the [link](https://progenomes.embl.de/data/proGenomes2.1_specI_clustering.tab) |
| proGenomes2.1_specI_lineageNCBI.tab | ProGenomes v.2 genome and taxonomy table, available in the resource through the [link](https://progenomes.embl.de/data/proGenomes2.1_specI_lineageNCBI.tab) |

