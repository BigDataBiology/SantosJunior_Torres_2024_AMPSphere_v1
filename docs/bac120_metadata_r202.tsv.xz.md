**bac120_metadata_r202.tsv**

**Description:**	GTDB v.202 Bacteria metadata. This TSV-like file describes the genomes in the
                        GTDB and gives their complete taxonomy as well as quality criteria.

It can be downloaded from the following [link][(https://data.gtdb.ecogenomic.org/releases/release202/202.0/)

| **column** | **description** | 
| :---: | :---: |
| accession | genome accession in the GTDB database |
| ambiguous_bases | number of 'Ns' in the genome |
| checkm_completeness | completeness in % |
| checkm_contamination | contamination in % |
| checkm_marker_count | number of marker genes identified |
| checkm_marker_lineage | lineage identified from the marker genes |
| checkm_marker_set_count | number of markers matching the lineage set |
| checkm_strain_heterogeneity | number of divergent markers belonging to the same lineage set |
| coding_bases | number of base pairs coding for proteins |
| coding_density | number of genes per base pair |
| contig_count | number of contigs in the genome |
| gc_count | number of G/C bases in the genome |
| gc_percentage | proportion of GC in the genome (%) |
| genome_size | length of the genome (base pairs) |
| gtdb_genome_representative | representative genome of the cluster in GTDB |
| gtdb_representative | is the representative? True/False |
| gtdb_taxonomy | full GTDB taxonomy lineage |
| gtdb_type_designation | is the type strain of species? True/False |
| gtdb_type_designation_sources | which is the type source? |
| gtdb_type_species_of_genus | is the type strain of the genus? True/False |
| l50_contigs | L50 of contigs |
| l50_scaffolds | L50 of scaffolds |
| longest_contig | length of longest contig |
| longest_scaffold | length of longest scaffolded sequence |
| lsu_23s_contig_len | length of contig containing the 23S rRNA gene |
| lsu_23s_count | number of 23S rRNA genes found |
| lsu_23s_length | length of 23S rRNA gene found |
| lsu_23s_query_id | access of the 23S rRNA gene found |
| lsu_5s_contig_len | length of contig containing the 5S rRNA gene |
| lsu_5s_count | number of 5S rRNA genes found |
| lsu_5s_length | length of 5S rRNA gene found |
| lsu_5s_query_id | access of the 23S rRNA gene found |
| lsu_silva_23s_blast_align_len | Blast alignment length of the search using 23S rRNA gene as query against SILVA database |
| lsu_silva_23s_blast_bitscore | Blast bitscore of the search using 23S rRNA gene as query against SILVA database |
| lsu_silva_23s_blast_evalue | Blast E-value of the search using 23S rRNA gene as query against SILVA database |
| lsu_silva_23s_blast_perc_identity | Blast identity (%) of the search using 23S rRNA gene as query against SILVA database |
| lsu_silva_23s_blast_subject_id | Blast best-hit of the search using 23S rRNA gene as query against SILVA database  |
| lsu_silva_23s_taxonomy |  Taxonomy of the best-hit of Blast search using 23S rRNA gene as query against SILVA database |
| mean_contig_length | average contig length in bp |
| mean_scaffold_length | average scaffold length in bp |
| mimag_high_quality | High-quality according to MiMAG criteria? True/False |
| mimag_low_quality | Low-quality according to MiMAG criteria? True/False |
| mimag_medium_quality | Medium-quality according to MiMAG criteria? True/False |
| n50_contigs | N50 of contigs |
| n50_scaffolds | N50 of scaffold |
| ncbi_assembly_level | Level of assembly as established by NCBI (Chromosome, plasmid) |
| ncbi_assembly_name | NCBI assembly access code |
| ncbi_assembly_type | Type of assembly as in NCBI (complete, partial) |
| ncbi_bioproject | BioProject in NCBI |
| ncbi_biosample | BioSample in NCBI |
| ncbi_contig_count | number of contigs |
| ncbi_contig_n50 | N50 of contigs |
| ncbi_country | Country of origin |
| ncbi_date | Date of deposition in the NCBI |
| ncbi_genbank_assembly_accession | Assembly accession data in GenBank |
| ncbi_genome_category | Supp. information about sources of this genome |
| ncbi_genome_representation | Representation of the genome in GenBank |
| ncbi_isolate | is an isolate microorganism? Yes/No |
| ncbi_isolation_source | Isolation source (material, environmental description) |
| ncbi_lat_lon | geographical coordinates to trace back the isolation/sample source |
| ncbi_molecule_count | number of molecules spotted |
| ncbi_ncrna_count | number of non-coding RNA detected |
| ncbi_organism_name | name of the organism according to NCBI taxid |
| ncbi_protein_count | number of predicted proteins from the genome |
| ncbi_refseq_category | is refseq or associated to any refseq? |
| ncbi_rrna_count | number of ribosomal RNAs detected |
| ncbi_scaffold_count | number of scaffolded sequences |
| ncbi_scaffold_l50 | L50 of scaffolded sequences |
| ncbi_scaffold_n50 | N50 of scaffolded sequences |
| ncbi_scaffold_n75 | N75 of scaffolded sequences |
| ncbi_scaffold_n90 | N90 of scaffolded sequences |
| ncbi_seq_rel_date | Date of sequence release in NCBI |
| ncbi_spanned_gaps | number of spanned gaps in the sequences |
| ncbi_species_taxid | Taxonomy identifier of the genome species in the NCBI database |
| ncbi_ssu_count | number of small ribosomal subunit genes in the genome |
| ncbi_strain_identifiers | identifiers for the strain of the genome |
| ncbi_submitter | organ/researcher submitting the sequence |
| ncbi_taxid | code related to the taxonomy group from NCBI |
| ncbi_taxonomy | full lineage according to NCBI taxonomy |
| ncbi_taxonomy_unfiltered | not controlled full lineage according to NCBI taxonomy |
| ncbi_total_gap_length | length of gaps introduced in the genome |
| ncbi_total_length | total genome length in the NCBI |
| ncbi_translation_table | translation table used in the genes prediction |
| ncbi_trna_count | number of transporter RNAs in the genome |
| ncbi_type_material_designation | is type material? |
| ncbi_ungapped_length | total length excluding gaps |
| ncbi_unspanned_gaps | total length excluding spanned gaps |
| ncbi_wgs_master | - |
| protein_count | number of predicted proteins |
| scaffold_count | number of scaffolded sequences |
| ssu_contig_len | length of the contig containing the predicted small ribosomal subunit genes |
| ssu_count | number of predicted small ribosomal subunit genes |
| ssu_gg_blast_align_len | Blast alignment length of the search using small ribosomal subunit gene as query against GreenGenes database |
| ssu_gg_blast_bitscore | Blast bitscore of the search using small ribosomal subunit gene as query against GreenGenes database |
| ssu_gg_blast_evalue | Blast E-value of the search using small ribosomal subunit gene as query against GreenGenes database |
| ssu_gg_blast_perc_identity | Blast identity (%) of the search using small ribosomal subunit gene as query against GreenGenes database |
| ssu_gg_blast_subject_id | Blast best-hit of the search using small ribosomal subunit gene as query against GreenGenes database |
| ssu_gg_taxonomy | taxonomy of the genome from the search of predicted small ribosomal subunit genes against GreenGenes database |
| ssu_length | length of the predicted small ribosomal subunit genes |
| ssu_query_id | access of the predicted small ribosomal subunit genes used to genome taxonomy identification |
| ssu_silva_blast_align_len | Blast alignment length of the search using small ribosomal subunit gene as query against SILVA database |
| ssu_silva_blast_bitscore | Blast bitscore of the search using small ribosomal subunit gene as query against SILVA database |
| ssu_silva_blast_evalue | Blast E-value of the search using small ribosomal subunit gene as query against SILVA database |
| ssu_silva_blast_perc_identity | Blast identity (%) of the search using small ribosomal subunit gene as query against SILVA database |
| ssu_silva_blast_subject_id | Blast best-hit of the search using small ribosomal subunit gene as query against SILVA database |
| ssu_silva_taxonomy | taxonomy of the genome from the search of predicted small ribosomal subunit genes against SILVA database |
| total_gap_length | total length of gaps included in the genome |
| trna_aa_count | number of different amino acids transported by the different tRNAs predicted from the genome |
| trna_count | number of transporter RNAs predicted in the genome |
| trna_selenocysteine_count | number of transporter RNAs predicted in the genome related to the selenocysteine anti-codon|


**MD5 SUM:**	7b893f37276dc04fb8a0912292ae8216

**Size (MBytes):**	40.063446044921875

**Content sample (first 5 items):**

accession	ambiguous_bases	checkm_completeness	checkm_contamination	checkm_marker_count	checkm_marker_lineage	checkm_marker_set_count	checkm_strain_heterogeneity	coding_bases	coding_density	contig_count	gc_count	gc_percentage	genome_size	gtdb_genome_representative	gtdb_representative	gtdb_taxonomy	gtdb_type_designation	gtdb_type_designation_sources	gtdb_type_species_of_genus	l50_contigs	l50_scaffolds	longest_contig	longest_scaffold	lsu_23s_contig_len	lsu_23s_count	lsu_23s_length	lsu_23s_query_id	lsu_5s_contig_len	lsu_5s_count	lsu_5s_length	lsu_5s_query_id	lsu_silva_23s_blast_align_len	lsu_silva_23s_blast_bitscore	lsu_silva_23s_blast_evalue	lsu_silva_23s_blast_perc_identity	lsu_silva_23s_blast_subject_id	lsu_silva_23s_taxonomy	mean_contig_length	mean_scaffold_length	mimag_high_quality	mimag_low_quality	mimag_medium_quality	n50_contigs	n50_scaffolds	ncbi_assembly_level	ncbi_assembly_name	ncbi_assembly_type	ncbi_bioproject	ncbi_biosample	ncbi_contig_count	ncbi_contig_n50	ncbi_country	ncbi_date	ncbi_genbank_assembly_accession	ncbi_genome_category	ncbi_genome_representation	ncbi_isolate	ncbi_isolation_source	ncbi_lat_lon	ncbi_molecule_count	ncbi_ncrna_count	ncbi_organism_name	ncbi_protein_count	ncbi_refseq_category	ncbi_rrna_count	ncbi_scaffold_count	ncbi_scaffold_l50	ncbi_scaffold_n50	ncbi_scaffold_n75	ncbi_scaffold_n90	ncbi_seq_rel_date	ncbi_spanned_gaps	ncbi_species_taxid	ncbi_ssu_count	ncbi_strain_identifiers	ncbi_submitter	ncbi_taxid	ncbi_taxonomy	ncbi_taxonomy_unfiltered	ncbi_total_gap_length	ncbi_total_length	ncbi_translation_table	ncbi_trna_count	ncbi_type_material_designation	ncbi_ungapped_length	ncbi_unspanned_gaps	ncbi_wgs_master	protein_count	scaffold_count	ssu_contig_len	ssu_count	ssu_gg_blast_align_len	ssu_gg_blast_bitscore	ssu_gg_blast_evalue	ssu_gg_blast_perc_identity	ssu_gg_blast_subject_id	ssu_gg_taxonomy	ssu_length	ssu_query_id	ssu_silva_blast_align_len	ssu_silva_blast_bitscore	ssu_silva_blast_evalue	ssu_silva_blast_perc_identity	ssu_silva_blast_subject_id	ssu_silva_taxonomy	total_gap_length	trna_aa_count	trna_count	trna_selenocysteine_count
GB_GCA_000013845.2	0	100	0	332	o__Clostridiales (UID1375)	124	0	2428396	82.0379664388356	4	835660	28.2309174592107	2960088	RS_GCF_000013285.1	f	d__Bacteria;p__Firmicutes_A;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium_P;s__Clostridium_P perfringens	not type material	none	f	1	1	2897393	2897393	2897393	10	2900	CP000312.1	2897393	9	296	CP000312.1-#9	2900	5356	0	100	CP000312.11694.14596	Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae;Clostridium;Clostridium perfringens SM101	740022	740022	t	f	f	2897393	2897393	Complete Genome	ASM1384v2	na	PRJNA12521	SAMN02604026	none	none	none	2006-07-26	GCA_000013845.2	none	full	none	none	none	4	0	Clostridium perfringens SM101	2620	na	0	none	none	none	none	none	2006/07/26	0	1502	0	SM101	TIGR	289380	d__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium;s__Clostridium perfringens	d__Bacteria;x__Terrabacteria group;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Clostridiaceae;g__Clostridium;s__Clostridium perfringens;x__Clostridium perfringens SM101	0	2921996	11	95	none	2921996	0	none	2659	4	2897393	10	none	none	none	none	none	none	1510	CP000312.1	1509	2787	0	100	CP000312.9992.11504	Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae;Clostridium sensu stricto 1;Clostridium perfringens SM101	0	20	95	1
GB_GCA_000016465.1	14	99.55	0.23	767	f__Pasteurellaceae (UID4932)	440	0	1584718	87.4070135513253	1	689682	38.0405279812291	1813033	RS_GCF_001457655.1	f	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Pasteurellaceae;g__Haemophilus;s__Haemophilus influenzae	not type material	none	f	1	1	1813033	1813033	1813033	6	3109	CP000671.1-#6	1813033	5	109	CP000671.1	3109	5729	0	99.936	CP031684.957063.960195	Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Pasteurellaceae;Haemophilus;Haemophilus influenzae	1813033	1813033	t	f	f	1813033	1813033	Complete Genome	ASM1646v1	na	PRJNA16400	SAMN02603084	none	none	none	2007-06-07	GCA_000016465.1	none	full	none	none	none	1	0	Haemophilus influenzae PittEE	none	na	27	none	none	none	none	none	2007/06/07	0	727	6	PittEE	Center for Genomic Sciences, Allegheny-Singer Research Institute	374930	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pasteurellales;f__Pasteurellaceae;g__Haemophilus;s__Haemophilus influenzae	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pasteurellales;f__Pasteurellaceae;g__Haemophilus;s__Haemophilus influenzae;x__Haemophilus influenzae PittEE	0	1813033	11	58	none	1813033	0	none	1777	1	1813033	6	none	none	none	none	none	none	1536	CP000671.1-#4	1536	2837	0	100	CP031684.831797.833346	Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Pasteurellaceae;Haemophilus;Haemophilus influenzae	0	20	58	1
GB_GCA_000024525.1	0	100	0.89	454	o__Cytophagales (UID2936)	336	25	7481814	88.1119617375894	9	4258276	50.148941417161	8491258	GB_GCA_000024525.1	t	d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Cytophagales;f__Spirosomaceae;g__Spirosoma;s__Spirosoma linguale	type strain of species	LPSN; DSMZ	t	1	1	8078757	8078757	8078757	4	2807	CP001769.1	none	0	none	none	2807	5184	0	100	CP001769.6728674.6731505	Bacteria;Bacteroidota;Bacteroidia;Cytophagales;Spirosomaceae;Spirosoma;Spirosoma linguale DSM 74	943473	943473	f	f	t	8078757	8078757	Complete Genome	ASM2452v1	na	PRJNA28817	SAMN00002598	none	none	none	2010-01-13	GCA_000024525.1	none	full	none	none	none	9	0	Spirosoma linguale DSM 74	6947	na	11	none	none	none	none	none	2010/01/13	0	108	4	DSM 74	US DOE Joint Genome Institute (JGI-PGF)	504472	d__Bacteria;p__Bacteroidetes;c__Cytophagia;o__Cytophagales;f__Cytophagaceae;g__Spirosoma;s__Spirosoma linguale	d__Bacteria;x__FCB group;x__Bacteroidetes/Chlorobi group;p__Bacteroidetes;c__Cytophagia;o__Cytophagales;f__Cytophagaceae;g__Spirosoma;s__Spirosoma linguale;x__Spirosoma linguale DSM 74	0	8491258	11	49	assembly from type material	8491258	0	none	7129	9	8078757	4	none	none	none	none	none	none	1503	CP001769.1	1495	2761	0	100	CP001769.6726840.6728334	Bacteria;Bacteroidota;Bacteroidia;Cytophagales;Spirosomaceae;Spirosoma;Spirosoma linguale DSM 74	0	20	47	0
GB_GCA_000152525.1	0	96.45	0.11	813	o__Pseudomonadales (UID4488)	308	0	5429460	87.2609346977394	124	4087884	66.502120221936	6222097	RS_GCF_001457615.1	f	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas;s__Pseudomonas aeruginosa	not type material	none	f	21	1	242903	6222097	6222097	1	226	CH482383.1	6222097	3	110	CH482383.1	226	411	2.33e-113	99.558	CP033432.3936150.3939058	Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas;Pseudomonas aeruginosa	49572	6222097	f	f	t	98844	6222097	Scaffold	ASM15252v1	na	PRJNA16170	SAMN02595210	124	98844	none	2006-01-11	GCA_000152525.1	none	full	none	none	none	0	0	Pseudomonas aeruginosa C3719	none	na	0	1	1	6222097	6222097	6222097	2006/01/11	123	287	0	C3719	Broad Institute	350704	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas;s__Pseudomonas aeruginosa	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas;x__Pseudomonas aeruginosa group;s__Pseudomonas aeruginosa;x__Pseudomonas aeruginosa C3719	75099	6222097	11	37	none	6146998	0	AAKV00000000.1	5785	1	none	0	none	none	none	none	none	none	none	none	none	none	none	none	none	none	75099	16	37	1
[...]
