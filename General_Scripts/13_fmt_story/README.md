# Effect of AMPs in the FMT outcome

FMT analysis brings the data for:

	(a) filtering of metagenome-assembled genomes;
	(b) amp screening with macrel;
	(c) count and collect AMP sequences per genome;
	(d) merge amp information with donor and recipient pre/post FMT;
	(e) fisher test using the number of MAGs with less than and greater
	    than x AMPs per genome for each outcome.

To reproduce the data and panels:

```
# to reproduce all features in the supplementary figure 1:
    $ python3 main.py
```

The scripts in *utils/* folder are explained below:
|**script** | **function** |
| :---: | :---: |
|*fmt_startup.py*| select quality controlled MAGs and submit them to Macrel | 
|*evaluation.py*| count and collect AMP sequences for each MAG and merge with metadata |
|*merger.py*| merge counts and AMP sequences from MAGs to their original samples per FMT event |
|*ftests.py*| count MAGs divided by outcome and AMP/MAG cutoffs, and submit the 2x2 matrix to Fisher's exact test |
|*plot_per_genome.py*| plot proportion of outcomes per each AMP/MAG in recipient pre/post FMT and donors | 


Input files were obtained from the study of
[Sebastian et al.](https://www.biorxiv.org/content/10.1101/2021.09.30.462010v1.supplementary-material)

Detailed inputs list:

| **Input file** | **Description** |
| :---: | :---: |
| dat.transmission.celio.Rdata | detailed information about the FMTs and their outcomes with the MAGs species per event |
| FMT_fastas.tar.gz | tar gzipped file containing the fasta files for the binned MAGs |
| FMT_META.metadata.all.v01.xlsx | table containing FMT metadata |
| FMT_META.metadata.per_FMT.v01.xlsx | table containing selected FMT metadata |
| mag.taxonomy.gtdb.raw.tsv | taxonomy complete lineages of MAGs |

