**ar122_metadata_r202.tsv**

**Description:**	GTDB v.202 Archea metadata. This TSV-like file describes the genomes in the
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

**MD5 SUM:**	d1ecd66de3a24a091a0a4a5306dd4f36

**Size (MBytes):**	0.6801910400390625

**Content sample (first 5 items):**

accession	ambiguous_bases	checkm_completeness	checkm_contamination	checkm_marker_count	checkm_marker_lineage	checkm_marker_set_count	checkm_strain_heterogeneity	coding_bases	coding_density	contig_count	gc_count	gc_percentage	genome_size	gtdb_genome_representative	gtdb_representative	gtdb_taxonomy	gtdb_type_designation	gtdb_type_designation_sources	gtdb_type_species_of_genus	l50_contigs	l50_scaffolds	longest_contig	longest_scaffold	lsu_23s_contig_len	lsu_23s_count	lsu_23s_length	lsu_23s_query_id	lsu_5s_contig_len	lsu_5s_count	lsu_5s_length	lsu_5s_query_id	lsu_silva_23s_blast_align_len	lsu_silva_23s_blast_bitscore	lsu_silva_23s_blast_evalue	lsu_silva_23s_blast_perc_identity	lsu_silva_23s_blast_subject_id	lsu_silva_23s_taxonomy	mean_contig_length	mean_scaffold_length	mimag_high_quality	mimag_low_quality	mimag_medium_quality	n50_contigs	n50_scaffolds	ncbi_assembly_level	ncbi_assembly_name	ncbi_assembly_type	ncbi_bioproject	ncbi_biosample	ncbi_contig_count	ncbi_contig_n50	ncbi_country	ncbi_date	ncbi_genbank_assembly_accession	ncbi_genome_category	ncbi_genome_representation	ncbi_isolate	ncbi_isolation_source	ncbi_lat_lon	ncbi_molecule_count	ncbi_ncrna_count	ncbi_organism_name	ncbi_protein_count	ncbi_refseq_category	ncbi_rrna_count	ncbi_scaffold_count	ncbi_scaffold_l50	ncbi_scaffold_n50	ncbi_scaffold_n75	ncbi_scaffold_n90	ncbi_seq_rel_date	ncbi_spanned_gaps	ncbi_species_taxid	ncbi_ssu_count	ncbi_strain_identifiers	ncbi_submitter	ncbi_taxid	ncbi_taxonomy	ncbi_taxonomy_unfiltered	ncbi_total_gap_length	ncbi_total_length	ncbi_translation_table	ncbi_trna_count	ncbi_type_material_designation	ncbi_ungapped_length	ncbi_unspanned_gaps	ncbi_wgs_master	protein_count	scaffold_count	ssu_contig_len	ssu_count	ssu_gg_blast_align_len	ssu_gg_blast_bitscore	ssu_gg_blast_evalue	ssu_gg_blast_perc_identity	ssu_gg_blast_subject_id	ssu_gg_taxonomy	ssu_length	ssu_query_id	ssu_silva_blast_align_len	ssu_silva_blast_bitscore	ssu_silva_blast_evalue	ssu_silva_blast_perc_identity	ssu_silva_blast_subject_id	ssu_silva_taxonomy	total_gap_length	trna_aa_count	trna_count	trna_selenocysteine_count
GB_GCA_000200715.1	2	99.03	0	145	k__Archaea (UID2)	103	0	1814397	88.7198386767109	6	1173335	57.3798631198789	2045086	GB_GCA_000200715.1	t	d__Archaea;p__Thermoproteota;c__Nitrososphaeria;o__Nitrososphaerales;f__Nitrosopumilaceae;g__Cenarchaeum;s__Cenarchaeum symbiosum	not type material	none	f	1	1	1281157	2045086	2045086	1	2966	DP000238.1	2045086	1	111	DP000238.1	2966	5478	0	100	DQ397549.30903.33889	Archaea;Crenarchaeota;Nitrososphaeria;Nitrosopumilales;Nitrosopumilaceae;Cenarchaeum symbiosum A	340809	2045086	t	f	f	1281157	2045086	Chromosome	ASM20071v1	na	PRJNA202	SAMN02744041	6	1281157	USA:Santa Barbara	2006-11-20	GCA_000200715.1	derived from environmental sample	full	none	10-20 meters below sea level off the coast of Santa Barbara; CA; USA	none	1	0	Cenarchaeum symbiosum A	2017	na	3	none	none	none	none	none	2006/11/20	5	46770	1	none	DOE Joint Genome Institute	414004	d__Archaea;p__Thaumarchaeota;c__;o__Cenarchaeales;f__Cenarchaeaceae;g__Cenarchaeum;s__Cenarchaeum symbiosum	d__Archaea;x__TACK group;p__Thaumarchaeota;o__Cenarchaeales;f__Cenarchaeaceae;g__Cenarchaeum;s__Cenarchaeum symbiosum;x__Cenarchaeum symbiosum A	229	2045086	11	45	none	2044857	0	none	1924	1	2045086	1	none	none	none	none	none	none	1470	DP000238.1	1470	2660	0	99.32	AF083072.21771.23243	Archaea;Crenarchaeota;Nitrososphaeria;Nitrosopumilales;Nitrosopumilaceae;Cenarchaeum;Cenarchaeum symbiosum	229	18	40	0
GB_GCA_000247545.1	0	100	0	174	k__Archaea (UID146)	136	0	2178950	88.8308628084079	2	1350357	55.0510004402916	2452920	GB_GCA_000247545.1	t	d__Archaea;p__Thermoproteota;c__Thermoproteia;o__Thermoproteales;f__Thermoproteaceae;g__Pyrobaculum;s__Pyrobaculum oguniense	type strain of species	LPSN; DSMZ	f	1	1	2436033	2436033	2436033	1	3014	CP003316.1	none	0	none	none	3014	5566	0	100	CP003316.1433240.1436302	Archaea;Crenarchaeota;Thermoprotei;Thermoproteales;Thermoproteaceae;Pyrobaculum;Pyrobaculum oguniense TE7	1226460	1226460	f	f	t	2436033	2436033	Complete Genome	ASM24754v1	na	PRJNA42375	SAMN02604076	none	none	Japan: Tsuetate Hot Spring; Kumamoto prefecture	2012-02-10	GCA_000247545.1	none	full	none	none	none	2	0	Pyrobaculum oguniense TE7	2835	na	3	none	none	none	none	none	2012/02/10	0	99007	0	TE7	UCSC-Lowelab	698757	d__Archaea;p__Crenarchaeota;c__Thermoprotei;o__Thermoproteales;f__Thermoproteaceae;g__Pyrobaculum;s__Pyrobaculum oguniense	d__Archaea;x__TACK group;p__Crenarchaeota;c__Thermoprotei;o__Thermoproteales;f__Thermoproteaceae;g__Pyrobaculum;s__Pyrobaculum oguniense;x__Pyrobaculum oguniense TE7	0	2452920	11	48	assembly from type material	2452920	0	none	2943	2	2436033	2	none	none	none	none	none	none	1167	CP003316.1	1166	2154	0	100	CP003316.1430990.1433182	Archaea;Crenarchaeota;Thermoprotei;Thermoproteales;Thermoproteaceae;Pyrobaculum;Pyrobaculum oguniense TE7	0	19	42	0
GB_GCA_000306725.1	0	99.84	0	228	p__Euryarchaeota (UID49)	153	0	2638662	85.872449246917	1	1370789	44.6108705210187	3072769	GB_GCA_000306725.1	t	d__Archaea;p__Halobacteriota;c__Methanosarcinia;o__Methanosarcinales;f__Methanosarcinaceae;g__Methanolobus;s__Methanolobus psychrophilus	type strain of species	LPSN	f	1	1	3072769	3072769	3072769	3	2914	CP003083.1	3072769	4	117	CP003083.1	2912	5378	0	100	CP003083.1449769.1452685	Archaea;Halobacterota;Methanosarcinia;Methanosarciniales;Methanosarcinaceae;Methanolobus;Methanolobus psychrophilus R15	3072769	3072769	t	f	f	3072769	3072769	Complete Genome	ASM30672v1	na	PRJNA74005	SAMN00739943	none	none	none	2012-10-18	GCA_000306725.1	none	full	none	none	none	1	0	Methanolobus psychrophilus R15	none	na	10	none	none	none	none	none	2012/10/18	0	420950	3	R15	Beijing Institute of Genomics	1094980	d__Archaea;p__Euryarchaeota;c__Methanomicrobia;o__Methanosarcinales;f__Methanosarcinaceae;g__Methanolobus;s__Methanolobus psychrophilus	d__Archaea;p__Euryarchaeota;x__Stenosarchaea group;c__Methanomicrobia;o__Methanosarcinales;f__Methanosarcinaceae;g__Methanolobus;s__Methanolobus psychrophilus;x__Methanolobus psychrophilus R15	0	3072769	11	52	none	3072769	0	none	3004	1	3072769	3	none	none	none	none	none	none	1473	CP003083.1	1472	2719	0	100	CP003083.1452985.1454458	Archaea;Halobacterota;Methanosarcinia;Methanosarciniales;Methanosarcinaceae;Methanolobus;Methanolobus psychrophilus R15	0	20	52	0
GB_GCA_000375685.1	51	73.65	2.21	174	k__Archaea (UID146)	136	66.67	1045795	86.8678419658906	176	385673	32.0368719789407	1203892	GB_GCA_000375685.1	t	d__Archaea;p__Thermoproteota;c__Thermoproteia;o__Gearchaeales;f__Gearchaeaceae;g__AAA261-N23;s__AAA261-N23 sp000375685	not type material	none	f	22	22	43826	43826	none	0	none	none	23413	1	99	AQVQ01000015.1	none	none	none	none	none	none	6840	6840	f	f	t	18548	18548	Contig	ASM37568v1	na	PRJNA165545	SAMN02256477	176	18548	none	2013-04-20	GCA_000375685.1	derived from single cell	full	none	none	none	0	none	crenarchaeote SCGC AAA261-N23	0	na	none	none	none	none	none	none	2013/04/20	0	1130305	none	SCGC AAA261-N23	DOE Joint Genome Institute	1130305	d__Archaea;p__Crenarchaeota;c__;o__;f__;g__;s__	d__Archaea;x__TACK group;p__Crenarchaeota;x__unclassified Crenarchaeota;s__crenarchaeote SCGC AAA261-N23	0	1203892	none	none	none	1203892	0	AQVQ00000000.1	1479	176	none	0	none	none	none	none	none	none	none	none	none	none	none	none	none	none	0	18	34	0
[...]
