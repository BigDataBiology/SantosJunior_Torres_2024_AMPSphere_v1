# Clonality and Quality

This analysis consists of:

	(a) quality description of AMP sets using different metrics;
	(b) the distribution of unique genes in the AMPSphere;
	(c) the length distribution of AMPs in AMPSphere;
	(d) the overlap percent among mammalian guts.

To reproduce the data and panels of figure S1:

```
# to reproduce all features in the supplementary figure 1:
    $ python3 main.py
```

The scripts in *utils/* folder are explained below:
|**script**|**outputs**|**figure panel**|
| :---: | :---: | :---: |
|*qualtest.py*| quality analysis involving all 5 tests | panel A |
|*enrich.py*| tests enrichment of quality-controlled AMPs in annotated sets | panel A |
|*ugenes.py*| counts the unigenes per AMP and plot a bar chart of it | panel B | 
|*AMPfeatures.py*| histogram of length of AMPs in AMPSphere | panel C |
|*hmammal_guts.py*| heatmap of AMP overlap content in % from different mammalian guts | panel D | 


Detailed inputs list:

| **Input file** | **Description** |
| :---: | :---: |
| AMPSphere_v.2021-03.fna.xz | fasta file contaning gene sequences for AMPs, currently allocated in Zenodo
| AMPSphere_v.2021-03.faa.gz | fasta file contaning peptide sequences, currently allocated in Zenodo
| gmsc_amp_genes_envohr_source.tsv | table containing the detailed information of genes| gmsc, amp, sample, source, specI, is_metagenomic, geographic_location, latitude, longitude, general envo name, environment_material |
| SPHERE_v.2021-03.levels_assessment.tsv.gz | table contaning peptide clusters association, currently allocated in Zenodo |
| antifam_100.tsv.gz | list of matches to Antifam from GMSC analysis |
| result.tsv.gz | list of genes spotted as potentially false in the coordinates analysis from GMSC |
| metaproteome_100.tsv.gz:  list of genes spotted as potentially true by appearing in metaproteomes from GMSC analysis |
| metaT_100.tsv.gz | list of genes spotted as potentially true by appearing in metatranscriptomes from GMSC analysis |
| rnacode_seq.tsv.xz | list of genes spotted as potentially true by RNACode from GMSC analysis |
| rnacode_seq_false.tsv.xz | list of genes spotted as potentially false by RNACode from GMSC analysis |
| dramp_candidates.txt | list of AMPs in AMPSphere with homologs in DRAMP v3 |
| annotated_candidates.txt | list of AMPs in AMPSphere with homologs to either DRAMP v3, starpepDB, SmProt, GMGC, STsORFs |

