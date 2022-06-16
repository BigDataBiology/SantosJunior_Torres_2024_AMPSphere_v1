# Suptable_s1.py

Internal script
Run:

```
python3 uscripts/suptable_s1.py
```

### : Description :

Create a supplementary table with information about the assemblies, gene
and AMP prediction. 

### : Inputs :

1. **data/samples-min500k-assembly-prodigal-stats.tsv.gz:**

Table from the resource used to build AMPSphere v2022-03 version.
It consists of the following columns:

    sample_accession - biosample
    study - study name
    study_accession - accession of the project in NCBI/ENA
    source - source of the sample
    human_filtering - filtering off human genome reads
    inserts_raw - total number of inserts
    inserts_hq - high-quality inserts
    inserts_filtered - number of inserts after filtering by reference
    assembly_total_length - total base pairs assembled
    assembly_number - number of contigs in the assembly
    assembly_mean_length - mean contig length
    assembly_longest - longest contig length in the assembly
    assembly_shortest - shortest contig length in the assembly
    assembly_N_count - counting Ns
    assembly_Gaps - gaps in assembly
    assembly_N50 - N50 of assembly discounting Ns
    assembly_N50n - N50 of assembly counting Ns
    assembly_N70 - N70 of assembly discounting Ns
    assembly_N70n - N70 of assembly counting Ns
    assembly_N90 - N90 of assembly discounting Ns
    assembly_N90n - N90 of assembly counting Ns
    prodigal_complete - number of ORFs with start or end codons
    prodigal_no_end - number of ORFs without end codon
    prodigal_no_start - number of ORFs without start codon
    prodigal_no_start_end - number of ORFs without start or end codons
    prodigal_smorfs - number of small ORFs predicted with Prodigal
    prodigal_total_orfs - total number of ORFs
    smORFs - number of complete small ORFs

'*' Table available in this [link](https://github.com/BigDataBiology/global_data/tree/master/freeze.v2)

2. **data/gmsc_amp_genes_envohr_source.tsv.gz:**

Table from AMPSphere v2022-03 version. It consists in a table linking the
GMSC10 genes and the AMPs in the final resource. This table also contains
general metadata information. Its columns are:
    
    gmsc - gene accession in GMSC10 resource
    amp - AMP accession in AMPSphere v2022-03 resource
    sample - biosample from which we could identify the gene
    source - contig taxonomy from GTDB or ProGenomes2
    specI - specI cluster
    is_metagenomic - True or False
    geographic_location - Country, Continent, Ocean, etc.
    latitude, longitude - geographical coordinates in decimals
    general_envo_name - environment classification
    environment_material - environment material, e.g.: stool

3. **data/reduced_metadata.tsv.gz**

Table generated during the metadata analysis.
It has the same data as the table `metadata.tsv.gz`.
Columns:

    sample_accession - biosample accession
    geographic_location - location as country, continent, ocean, etc.
    latitude - decimal latitude coordinate
    longitude - decimal longitude coordinate
    general_envo_name - habitat
    environment_material - sampled material, e.g. stool

### : Outputs :

1. **supplementary_table_S1.tsv**

Table with supplementary information by environment.
Columns:

    sample - biosample accession
    habitat - environment
    raw inserts - number of inserts in a sample
    assembled bp - assembled base pairs
    N50 - N50 of assembly
    ORFs+smORFs - total number of ORFs in a given environment, counting small ORFs
    smORFs - number of small ORFs in a given environment
    non-redundant AMPs - number of non-redundant AMPs in a given environment
    
