# Scripts needed to reproduce Figure 1 of AMPSphere paper:

AMPSphere collates 836,498 non-redundant AMP candidates from thousands of
metagenomes and high-quality microbial genomes. It consists of:

	(a) the geographical distribution of the publicly available metagenomes
	from diverse sources;

	(b) the limited portion of AMP candidates with homologs in peptide
	(SmProt, DRAMP, starPepDB, STsORFs) and large protein datasets - GMGC;

	(c) a small overlap of AMPs belonging to different environments;

	(d) overlap of AMPs among human body sites. 

The code can be executed as follows:

```
# to reproduce all features in the figure 1:
    $ python3 main.py
```

The scripts in *utils/* folder are explained below:
|**script**|**inputs**|**outputs**|**figure panel**|
| :---: | :---: | :---: | :---: |
|*download_files.py*| - | folder _data_ contaning all files needed | - |
|*mundi_map.py*| gmsc_amp_genes_envohr_source.tsv | map of samples distribution containing AMPs | panel A |
|*homologs.py*| AMPSphere.faa from Zenodo repository and all files in data/databases_homology | generates several info | panel B |
|*henvo.py*| gmsc_amp_genes_envohr_source.tsv | heatmap of AMP overlap content in % from different environments | panel C |
|*hbody.py*| gmsc_amp_genes_envohr_source.tsv | heatmap of AMP overlap content in % from different human body sites | panel D | 

Detailed inputs list:

 - AMPSphere_v.2021-03.faa.gz:	fasta file contaning peptide sequences, currently allocated in Zenodo

 - gmsc_amp_genes_envohr_source.tsv:	table containing the detailed information of genes - gmsc, amp,
 sample, source, specI, is_metagenomic, geographic_location, latitude, longitude, general envo name,
 environment_material

 - all_SmProt.fa.gz:	[Small proteins database (SmProt v.2)](http://bigdata.ibp.ac.cn/SmProt/)

 - starPepDB.fasta:	[starPep database](http://mobiosd-hub.com/starpep/)

 - DRAMP.fa:	[DRAMP v.3](DRAMP.cpu-bioinfor.org/)

 - STsORFs.faa: [Salmonella typhimurium small proteins catalog](https://www.biorxiv.org/content/10.1101/2020.05.26.116038v1)

 - GMGC10.proGenomes.faa:	[Global microbial genes catalog v.1](gmgc.embl.de/)

 - truepep_gmgc_progenomes.m8.xz:	pre-computed table of hits searching AMPSphere against GMGC
 
