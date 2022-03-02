# Figure 2 - AMP candidates are mostly accessory from a small number of microbial species

Figure 2 brings data about the:

	(a) taxonomic affiliation of AMPs in our database;
	(b) the distribution of AMPs within the Prevotella genus - the more abundant source of AMPs in the AMPSphere;
	(c) quality families span a small number of genera;
	(d) AMP candidates and families from ProGenomes2 along with their high quality subsets were classified as core, shell, and accessory genes.

To reproduce the data and panels of figure 2:

```
# to reproduce all features in the supplementary figure 1:
    $ python3 main.py
```

The scripts in *utils/* folder are explained below:
|**script**|**outputs**|**figure panel**|
| :---: | :---: | :---: |
|*download_files.py*| folder _data_ contaning all files needed | - |
|*taxonomy.py*| graph of taxonomy affiliation of AMPs in AMPSphere | panel A |
|*prevotella_amp.py*| a table of AMPs within *Prevotella* genus, and a tree of *Prevotella* species | panel B |
|*family_sizes.py*| a graph of the number of genera per quality-controlled family | panel C | 
|*core.py*| bar chart of the amount of genes associated as core, shell and strain-specific | panel D |

Detailed inputs list:

 - gmsc_amp_genes_envohr_source.tsv:	table containing the detailed information of genes - gmsc, amp, sample,
 source, specI, is_metagenomic, geographic_location, latitude, longitude, general envo name,
 environment_material

 - SPHERE_v.2021-03.levels_assessment.tsv.gz:	table contaning peptide clusters association, currently allocated in Zenodo

 - proGenomes2.1_specI_clustering.tab:    ProGenomes v.2 specI and samples table, available in the resource through
 the [link](https://progenomes.embl.de/data/proGenomes2.1_specI_clustering.tab)
 
 - quality_families.txt:    table of families containing: family, has_experimental_evidence,
 number_of_quality_candidates, total_AMPs, percent of quality-controlled candidates

 
 - quality_candidates.txt:    list of AMPs passing all tests without experimental evidence of transcription/translation
 
 - high_quality_candidates.txt:    list of AMPs passing all tests and also presenting experimental evidence of
 transcription/translation
 
 - prevotella_species_list.tsv:    curated table of *Prevotella* species with habitats, hosts and info enough
 to build the iToL metadata from [LPSN](https://lpsn.dsmz.de/)
 
 - 37185.fasta:    Fasta file containing the 16S sequences for the *Prevotella* species,
 curated by [LPSN](https://lpsn.dsmz.de/)
 
 - complete_amps_associated_taxonomy.tsv.gz:    complete table of annotation of AMPs and their genes containing:
 gmsc, amp, sample, contig, start, stop, strand, ID, partial, start_type, rbs_motif, rbs_spacer, gc_cont, taxid,
 level, name, retained, assigned, agreement, support, specI
 
