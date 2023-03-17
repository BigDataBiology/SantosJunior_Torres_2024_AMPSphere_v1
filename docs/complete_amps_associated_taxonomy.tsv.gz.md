**complete_amps_associated_taxonomy.tsv**

**Description:**	TSV-file containing the metagenome candidate AMP genes description, with the information obtained
                        from their prediction using Prodigal_sm.

| **columns** | **description** |
| :---: | :---: |
| amp | AMPSphere candidate AMP access code |
| gmsc | candidate AMP gene access code in the GMSC resource |
| sample | metagenome sample of origin |
| contig | contig name |
| start | small open reading frame start position in the contig (nucleotide) |
| stop | small open reading frame stop position in the contig (nucleotide) |
| strand | strand sense (+1 or -1) |
| fid | random gene identifier associated inside the gene prediction for each metagenome |
| partial | is partial? (contains start and stop - 0, contains only start or stop - 1)|
| start_type | start codon |
| rbs_motif | motif associated to the ribosome binding site |
| rbs_spacer | length of the spacer left between the ribosome binding site and the gene start codon |
| gc_cont | content of GC in the small open reading frame (%) |
| taxid | taxonomy group using TAXID from NCBI after contig identification using MMSeqs2 approach with GTDB database |
| level | lowest taxonomy rank achieved in the contig identification using MMSeqs2 approach with GTDB database |
| source | taxonomy group name/microbial candidate AMP source |
| retained | - |
| assigned | - |
| agreement | - |
| support | confidence in the taxonomy identification (%) |
| specI | species cluster according to ProGenomes v2 |

**MD5 SUM:**	8762d63a7acceff48c84ad644241c99f

**Size (MBytes):**	183.7403974533081

**Content sample (first 5 items):**

amp	gmsc	sample	contig	start	stop	strand	fid	partial	start_type	rbs_motif	rbs_spacer	gc_cont	taxid	level	source	retained	assigned	agreement	support	specI
AMP10.000_000	GMSC10.SMORF.000_036_844_556	SAMEA104142073	k119_56222	123	209	1	82432_1	0	ATG	TAA	7bp	0.333	237.0	genus	Phocaeicola	70.0	61.0	39.0	0.76	*
AMP10.000_000	GMSC10.SMORF.000_036_899_109	SAMEA104142074	k119_94475	12282	12368	1	76658_13	0	ATG	TAA	7bp	0.333	238.0	species	Phocaeicola vulgatus	36.0	36.0	12.0	0.53	specI_v3_Cluster2367
AMP10.000_000	GMSC10.SMORF.000_036_944_928	SAMEA104142075	k119_139920	12248	12334	1	89279_13	0	ATG	TAA	7bp	0.333	238.0	species	Phocaeicola vulgatus	39.0	39.0	12.0	0.53	specI_v3_Cluster2367
AMP10.000_000	GMSC10.SMORF.000_036_987_159	SAMEA104142076	k119_139643	10084	10170	-1	65412_13	0	ATG	TAA	7bp	0.322	237.0	genus	Phocaeicola	23.0	23.0	10.0	0.68	*
[...]
