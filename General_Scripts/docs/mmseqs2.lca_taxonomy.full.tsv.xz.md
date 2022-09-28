**mmseqs2.lca_taxonomy.full.tsv**

**Description:**	TSV-file with the taxonomy of contigs assembled for the metagenomes used in the
                        AMPSphere resource. The taxonomy was identified using the [MMSeqs2 procedure](https://github.com/soedinglab/mmseqs2/wiki#taxonomy-assignment)
                        with the GTDB database r.95 as reference.

| **columns** | **description** |
| :---: | :---: |
| sample | metagenome sample |
| contig | contig name |
| taxid | taxonomy identifier as in the NCBI |
| level | taxon rank (kingdom, genus, etc.) |
| name | taxon name (binomial when available) |
| retained | number of markers found significant |
| assigned | number of markers assigned to any taxonomy |
| agreement | number of markers agreeing to the taxonomy assigned |
| support | statistical support for the taxonomy assigned to the contig |


**MD5 SUM:**	de00608ab2e13cc60bb8060b7617ea6d

**Size (MBytes):**	59.12626266479492

**Content sample (first 5 items):**

sample	contig	taxid	level	name	retained	assigned	agreement	support
Karasov_2018_arabidopsis_NextMet1	k141_9474	295	genus	Sphingomonas	3	3	2	0.66
Karasov_2018_arabidopsis_NextMet124	k141_35350	31575	species	Microcoleus sp000317475	1	1	1	1.0
Karasov_2018_arabidopsis_NextMet124	k141_21720	33961	species	Sphingomonas sp002796605	1	1	1	1.0
Karasov_2018_arabidopsis_NextMet124	k141_42331	69	species	Pseudomonas E viridiflava	4	4	3	0.98
[...]
