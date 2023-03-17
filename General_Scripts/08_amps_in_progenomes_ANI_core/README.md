## AMPs in progenomes: ANI and Core analysis

This scripts set will download the genomes and proteomes of isolates from each
one of the species clusters in ProGenomes v2. After that, it will calculate the
ANI - Average Nucleotide Identity for each pair of genomes within the clusters, 
will cluster them as clones (ANI > 99.99%) and will separate them into
strains (99.5 < ANI < 99.99). The proteomes will be stored as non-redundant proteins, 
that will be reduced using a shortened alphabet and later on clustered using CDHIT
similar to the clustering procedures adopted in AMPSphere. Then, the prevalence of
protein families of at least 8 protein variants will be evaluated.

The major instructions to run the code in this folder are in the `main.sh` script.
of this analysis, the user should simply:

```
    $ ./main.sh
```

The inputs used are all described bellow:

| **Input file** | **Description** |
| :---: | :---: |
| proGenomes2.1_specI_lineageNCBI.tab.xz | table from ProGenomes2.1, with genome sample, kingdom, phylum, class, order, family, genus, species according NCBI standards |
| proGenomes2.1_specI_clustering.tab.xz | table from ProGenomes2.1, with species cluster and genome sample | 
| SPHERE_v.2022-03.levels_assessment.tsv.gz | table generated in the AMPSphere resource, it consists of the AMPs and their clustering information |
| complete_gmsc_pgenomes_metag.tsv.gz | table generated during metadata analysis with columns amp, gene, sample, source (species), is_metagenomic (bool: True/False) |

The outputs generated consist of:

| **Output file** | **Description** |
| :---: | :---: |
| amp_results_ani.tsv.xz | output table with the columns species, genome1, genome2, strain1, strain2, shared_amps, total_nr_amps, percent_shared -- Used in the calculation of odds ratio for pairs of genomes from the same strain to share AMPs |
| summary_output_core_prots.tsv.xz | output table with the columns cluster, core_fams, shell_fams, accessory_fams -- It brings the percent of full-length protein families calculated as in the AMPSphere belonging to each class for each species cluster |

