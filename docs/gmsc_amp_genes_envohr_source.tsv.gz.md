**gmsc_amp_genes_envohr_source.tsv**

**Description:**	Summarized information of candidate AMPs found in the metagenomes and progenomes.
                        This TSV-file contains most of the metadata used througout AMPSphere analysis.

| **columns** | **description** |
| :---: | :---: |
| gmsc | small ORF unique identifier in the GMSC resource |
| amp | candidate AMP in the AMPSphere resource |
| sample | Metagenome/ProGenomes v2 sample |
| source | Taxonomic group of the contig from which the AMP was identified |
| specI | species cluster as defined in ProGenomes v2 |
| is_metagenomic | if originated from metagenome - True, else condition - False |
| geographic_location | location from which the sample was taken (name) |
| latitude | decimal geographic coordinate of latitude |
| longitude | decimal geographic coordinate of longitude |
| general_envo_name | high-level habitat group |
| environmental_material | material used to isolate the sample, e.g. stool | 

NOTE: candidate AMPs from ProGenomes v2 samples were not included in the metadata, so for this table, they do not present the following fields:
      geographic_location, latitude, longitude, general_envo_name, environment_material -- These fields were left as 'NA'

**MD5 SUM:**	4447ee32a0a0db76da35f79df1daecdb

**Size (MBytes):**	110.56013774871826

**Content sample (first 5 items):**

gmsc	amp	sample	source	specI	is_metagenomic	geographic_location	latitude	longitude	general_envo_name	environment_material
GMSC10.SMORF.000_000_962_653	AMP10.000_000	435590.SAMN02604309	Bacteroides vulgatus	specI_v3_Cluster2367	False	NA	NA	NA	NA	NA
GMSC10.SMORF.000_001_397_691	AMP10.000_000	457394.SAMN02463694	Bacteroides sp. 4_3_47FAA	specI_v3_Cluster2367	False	NA	NA	NA	NA	NA
GMSC10.SMORF.000_001_600_202	AMP10.000_000	469593.SAMN02463762	Bacteroides sp. 3_1_40A	specI_v3_Cluster2367	False	NA	NA	NA	NA	NA
GMSC10.SMORF.000_001_968_479	AMP10.000_000	997891.SAMN02463934	Bacteroides vulgatus	specI_v3_Cluster2367	False	NA	NA	NA	NA	NA
[...]
