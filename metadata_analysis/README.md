# Metadata analysis

This analysis integrates all resources from AMPSphere and is a crucial part of the web platform.
It intends to generate a set of tables that besides useful for analysis, including some in the
figures that made up to the paper, also are used in the data management.

The set of analysis include:

	(a) creating of a convertion mechanism to get specI clusters from GTDB taxonomic affiliation;
	(b) get organized tables of genes from GMSC with their info, such as contigs, samples, start and stop positions, gc_content and etc.
	(c) add the taxonomy annotation using GTDB to the genes table
	(d) normalize the genes info adding specI clusters and AMP association
	(e) add metadata info

To reproduce the analysis:

```
# to reproduce all features in the supplementary figure 1:
    $ python3 main.py
```

It generates three main tables:

    I. pgenomes_AMP_specI.tsv:
    
    This table contains the AMP genes predicted from genomes in ProGenomes v2.
    The columns are the following:
    
    |**Column**|**Description**|
    | :---: | :---: |
    | AMP | Antimicrobial peptide accession in AMPSphere |
    | GMSC10 | Gene accession in GMSC resource |
    | genome | genome accession in ProGenomes v2 database |
    | start | gene start position in the contig/genome |
    | stop | gene stop position in the contig/genome |
    | species | NCBI taxonomy (as in ProGenomes resource) |
    | specI | specI cluster (genomic species) |

    II. complete_gmsc_pgenomes_metag.tsv:
    
    This table is a full resource of genes and their basic info, containing their origin and their samples.
    The columns are as follow:

    |**Column**|**Description**|
    | :---: | :---: |
    | gmsc | Gene accession in GMSC resource |
    | amp | Antimicrobial peptide accession in AMPSphere |
    | sample | biosample accession number |
    | source | taxonomy affiliation |
    | specI | specI cluster (genomic species) |
    | is_metagenomic | True (comes from metagenomes), False (comes from ProGenomes v2) |

    III. gmsc_amp_genes_envohr_source.tsv:
    
    This table is the full resource produced by the present analysis. It was used in the AMPSphere website.
    It contains the following columns:

    |**Column**|**Description**|
    | :---: | :---: |
    | gmsc | Gene accession in GMSC resource |
    | amp | Antimicrobial peptide accession in AMPSphere |
    | sample | biosample accession number |
    | source | taxonomy affiliation |
    | specI | specI cluster (genomic species) |
    | is_metagenomic | True (comes from metagenomes), False (comes from ProGenomes v2) |
    | geographic_location | Country/Geographical area from where sample comes from |
    | latitude | Latitude in decimals |
    | longitude | Longitude in decimals |
    | general envo name | Human readable environment name |
    | environment_material | Material sampled |
    
The scripts in *utils/* folder are explained below:

	- download_files.py:     **Function.** downloads the set of inputs needed for the other scripts
	                         **Inputs.** No inputs are needed 
	                         **Outputs.** Folder *data/* with diverse needed files (more details about the files below)

	- gtdb2pgenomes.py:      **Function.** Generates a conversion mechanism from GTDB taxonomy into specI clusters based in the samples
	                         **Inputs.** Resources from GTDB (bac120_metadata_r202.tsv, ar122_metadata_r202.tsv) and from ProGenomes (progenomes_samples.tsv).
	                         **Outputs.** Two tables: **(i)** gtdb_to_pgenomes.tsv (metadata of genomes in GTDB that can be converted into specI: "accession",
	                         "checkm_marker_lineage", "gtdb_taxonomy", "ncbi_assembly_name", "ncbi_bioproject", "ncbi_biosample", "ncbi_organism_name",
	                         "ncbi_taxid", "ncbi_taxonomy", "cluster"), and **(ii)** conv_pgen_to_gtdb_list.tsv (a complete association table of sample
	                         and the lineage - columns: "ncbi biosample", "specI cluster", "Domain", "Phylum", "Class", "Order", "Family", "Genus",
	                         and "Species")

	- metadata.py:           **Function.** Reduces the metadata generated during sampling and convert the name of the environments to human-readable format
	                         **Inputs.** metadata.tsv and general_envo_names.tsv
	                         **Outputs.** Generates the table reduced_metadata.tsv (columns: "sample_accession", "geographic_location", "latitude",
	                         "longitude", "general envo name", "environment_material")

	- progenomes_genes.py:   **Function.** Generates a table of AMP genes from ProGenomes by using their coordinates of gene prediction.
	                         **Inputs.** AMP genes list from AMPSphere_origins in Zenodo repository, the ampdictionary made from the genes fasta file in
	                         the Zenodo repository and the pre-computed resource from GMSC (GMSC10.ProGenomes2.coords.txt.gz)
	                         **Outputs.** Table AMPSphere_ProGenomes2.tsv.gz (columns: AMP, GMSC10, genome, contig, start, stop)

	- preprocess.py:         **Function.** Pre-process inputs from taxonomic annotation using GTDB, gene predicition and AMP classification. 
	                         **Inputs.** Pre-computed for GTDB contigs annotation (mmseqs2.lca_taxonomy.full.tsv.xz),
	                         gene prediction in GMSC (GMSC10.metag_smorfs.rename.txt.xz),
	                         and AMPSphere files from Zenodo (AMPSphere_v.2021-03.fna.xz,  AMPSphere_v.2021-03.origin_samples.tsv.gz)
	                         **Outputs.** Tables: GMGCgenes, gmsc_genes_smorfs.txt, gmsc_metag_amps.tsv.gz,
	                         amp_contigs_filtered_mmseqs2.lca_taxonomy.tsv.xz, gmsc_meta_taxo.tsv.gz

	- addspecI.py:           **Function.** Add specI clusters into complete gene tables.
	                         **Inputs.** gmsc_meta_taxo.tsv.gz and conv_pgen_to_gtdb_list.tsv
	                         **Outputs.** Table complete_amps_associated_taxonomy.tsv.gz (columns: gmsc, amp, sample, contig, start, stop, strand, ID,
	                         partial, start_type, rbs_motif, rbs_spacer, gc_cont, taxid, level, name, retainedassigned, agreement, support, specI)

	- hr_envo.py:            **Function.** Add metadata to the pre-computed gene table and concatenate the genomes from metagenomes and progenomes.
	                         **Inputs.** pre-computed files (AMPSphere_ProGenomes2.tsv.gz, complete_amps_associated_taxonomy.tsv.gz, reduced_metadata.tsv),
	                         and progenomes resource (pgenomes_samples.tsv)
	                         **Outputs.** It generates the three main tables previously described


Detailed inputs list:

 - proGenomes2.1_specI_clustering.tab:    ProGenomes v.2 specI and samples table, available in the resource through
 the [link](https://progenomes.embl.de/data/proGenomes2.1_specI_clustering.tab)

 - bac120_metadata_r202:    GTDB bacterial genomes metadata, available in the resource through
 the [link](https://data.gtdb.ecogenomic.org/releases/release202/202.0/bac120_metadata_r202.tar.gz)

 - ar122_metadata_r202:    GTDB archeal genomes metadata, available in the resource through
 the [link](http://data.gtdb.ecogenomic.org/releases/release202/202.0/ar122_metadata_r202.tar.gz)

 - AMPSphere_v.2021-03.fna.xz:	fasta file contaning gene sequences for AMPs, currently allocated in Zenodo

 - AMPSphere_v.2021-03.origin_samples.tsv.gz:	Table containing the AMP, genes from GMSC associated with it,
 and the metagenome samples and genomes in which it was spotted,
 currently allocated in Zenodo

 - metadata.tsv:	Table containing the metadata for the metagenome samples used to predict smORFs, currently
 it is allocated in the [public repository](https://github.com/BigDataBiology/global_data/tree/master/freeze.v2)

 - general_envo_names.tsv:	Table of associations of the original environment name and the general and humand readable
 ones. It is allocated in the [public repository](https://github.com/BigDataBiology/global_data/tree/master/freeze.v2)

 - GMSC10.ProGenomes2.coords.txt.gz:	Internal GMSC resource. Table with coordinates for all predicted smORFs in the
 genomes from ProGenomes2

 - GMSC10.metag_smorfs.rename.txt.xz:	Internal GMSC resource. Table with all coordinates and predictions of GMSC genes
 from metagenomes, with accession codes to GMSC database

 - mmseqs2.lca_taxonomy.full.tsv.xz:	Internal GMSC resource. Annotation of contigs taxonomy signature by using mmseqs2
 and the GTDB database

