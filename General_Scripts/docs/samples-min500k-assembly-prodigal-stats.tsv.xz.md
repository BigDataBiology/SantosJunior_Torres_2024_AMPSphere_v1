**samples-min500k-assembly-prodigal-stats.tsv**

**Description:**	Supplemmentary metadata for the assemblies used in the AMPSphere. 
                        This TSV-file consists of a set of statistics for each metagenome assembled sample.

| **columns** | **description** |
| :---: | :---: |
| sample_accession | metagenome sample accession in biosample db |
| study | study name |
| study_accession | bioproject accession on the NCBI |
| source | database from which the sequences come from |
| human_filtering | filtered against host genome? (human, or other species?) |
| inserts_raw | number of reads |
| inserts_hq | number of high-quality reads |
| inserts_filtered | number of remaining reads after quality-control |
| assembly_total_length | total length of the assembly in bp |
| assembly_number | number of contigs |
| assembly_mean_length | mean contig length |
| assembly_longest | length of the longest contig |
| assembly_shortest | length of the shortest contig |
| assembly_N_count | number of N's added to the assembly |
| assembly_Gaps | number of gaps in the contigs |
| assembly_N50 | N50 excluding Ns |
| assembly_N50n | N50 assuming Ns |
| assembly_N70 | N70 excluding Ns |
| assembly_N70n | N70 assuming Ns |
| assembly_N90 | N90 excluding Ns |
| assembly_N90n | N90 assuming Ns |
| prodigal_complete | number of complete ORFs from Prodigal prediction |
| prodigal_no_end | number of ORFs missing stop codon from Prodigal prediction |
| prodigal_no_start | number of ORFs missing start codon from Prodigal prediction |
| prodigal_no_start_end | number of ORFs without start and stop codons from Prodigal prediction |
| prodigal_smorfs | Number of ORFs predicted by prodigal that are shorter than 90 nucleotides. Note that these are not included in the above counts as they are not used in GMGC building |
| prodigal_total_orfs | total number of ORFs predicted from the assembly regardless completeness or size |
| smORFs | Number of ORFs predicted by the modified version of prodigal shipped with macrel. Only complete ORFs are predicted. This is a separate process from the one that generates prodigal_smorfs. In fact, ORFs counted in prodigal_smorfs are always incomplete: prodigal will generate complete ORFs with a minimum size of 90, but incomplete ones starting at 60 |

**MD5 SUM:**	c6b090f2b1062a409eec8c1cb3f994dd

**Size (MBytes):**	3.6716766357421875

**Content sample (first 5 items):**

sample_accession	study	study_accession	source	human_filtering	inserts_raw	inserts_hq	inserts_filtered	assembly_total_length	assembly_number	assembly_mean_length	assembly_longest	assembly_shortest	assembly_N_count	assembly_Gaps	assembly_N50	assembly_N50n	assembly_N70	assembly_N70n	assembly_N90	assembly_N90n	prodigal_complete	prodigal_no_end	prodigal_no_start	prodigal_no_start_end	prodigal_smorfs	prodigal_total_orfs	smORFs
Karasov_2018_arabidopsis_NextMet1	Karasov_2018_arabidopsis	no_id	LOCAL	Arabidopsis	12374844	12280133	2816546	6714448	10009	670.84	28303	201	0	0	700	2304	484	4634	348	7932	1908	2947	2852	3276	565	10983	2254
Karasov_2018_arabidopsis_NextMet124	Karasov_2018_arabidopsis	no_id	LOCAL	Arabidopsis	17119714	16975595	4970623	29141743	43811	665.17	18651	200	0	0	684	10111	479	20408	348	34791	7674	14216	14227	17616	2086	53733	10438
Karasov_2018_arabidopsis_NextMet15	Karasov_2018_arabidopsis	no_id	LOCAL	Arabidopsis	13384712	13259192	3511093	3583114	3780	947.91	19886	247	0	0	1325	644	717	1387	401	2740	1565	933	895	267	328	3660	1046
Karasov_2018_arabidopsis_NextMet25	Karasov_2018_arabidopsis	no_id	LOCAL	Arabidopsis	15623952	15495611	2985496	5459772	7847	695.78	16929	202	0	0	736	1622	486	3486	350	6153	1670	1986	2022	2116	514	7794	1998
[...]
