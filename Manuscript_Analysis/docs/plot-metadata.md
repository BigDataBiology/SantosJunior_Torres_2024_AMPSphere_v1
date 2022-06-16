# Plot-metadata.py

Internal script
Run:

```
python3 uscripts/plot-metadata.py
```

### : Description :

Plots and test the differences of AMP density (non-redundant AMPs per assembled
Mbps) in host-associated samples.

### : Inputs :

1. **data/gmsc_amp_genes_envohr_source.tsv.gz:**

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

2. **data/samples-min500k-assembly-prodigal-stats.tsv.gz:**

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

3. **data/metadata.tsv.gz:**

Table containing metadata description of the samples used in AMPSphere v2022-03 version.
This table can be found in the [link](https://github.com/BigDataBiology/global_data/tree/master/freeze.v2).
It consists of the following columns:

    sample_accession - biosample accession
    ena_ers_sample_id - name of the sample
    database - database available 
    access_status - status as 'public' or not
    study - study name
    study_accession - accession of study 
    publications - reference
    collection_date - date
    aliases - alternative accession
    microontology - microontology classification
    environment_biome - biome ontology
    environment_featureenvironment_material - sampled material, e.g. stool
    geographic_location - location of the sample: country, continent, ocean...
    latitude - decimal latitude
    longitude - decimal longitude
    tax_id - NCBI taxID of organism
    tax_scientific_name - taxonomy of organism
    host_common_name - host common name
    host_scientific_name - host scientific name
    host_tax_id - NCBI taxID of host
    gender - host gender
    age_years - host age
    subject_disease_status - details about disease
    antibiotic - use of antibiotics, details
    notes - commentary

4. **data/general_envo_names.tsv.gz:**

Table containing conversion tools to transform the microontology into
72 environments general habitats:

    # samples - number of samples
    # amps - number of amps detected in that habitat
    microontology - microontology term
    host_scientific_name - host species name
    host_tax_id - host NCBI taxID
    general_envo_name - general habitat name

### : Outputs :

1. **map_habitat.png**

Map containing the distribution of samples (as colored dots) per environment.

2. **samples_per_habitat.svg**

Bar chart containing the number of samples per habitat.

3. **hq_insters_per_habitat.svg**

Bar chart containing the number of inserts per habitat.

4. **smORFs_raw_millions.svg**

Raw number of smORFs per habitat as a horizontal bar chart.

5. **host_vs_nonhost.svg**

Box plot of host-associated and non-host associated samples in the `x axis`,
with the AMP density (number of non-redundant AMPs per sample normalized by
the assembled megabase pairs). Sampled dots were also plotted.

6. **ampsphere_amps_per_assembly_mbps.{svg, png}**

Box plot with `x axis` as the habitats, and the `y axis` as the AMP density as in 5.
Environments must contain at least 100 samples and AMPs. Habitats were colored as
host- and non-host-associated too.

7. **smorfs_per_assembly_mbps.svg**

Histogram of the number of small ORFs per assembled megabase pairs.

8. **amp_smorfs_sample.svg**

KDE plot with two curves:

    (1) number of small ORFs normalized by the assembled megabase pairs
    (2) number of non-redundant AMPs in AMPSphere normalized by
        the assembled megabase pairs

9. **mannwhitneyu_test_mammalguts.tsv**

Table reporting the Mann-Whitney U test performed pairwise in the
mammalian samples from guts. It was tested the number of non-redundant
AMPs per assembled megabase pairs per sample.

The mammalian species chosen need to have at 
least 100 samples and AMPs. It has the following columns:

    s1 - mammalian species 1
    s2 - mammalian species 2
    U - test statistics
    pval - probability, p-value
    log(pval) - log10 of p-value

10. **table_S2.tsv**

Table contaning the samples used in the AMPSphere study.
Columns:

    sample_accession - biosample accession
    ena_ers_sample_id - biosample accession
    database - database where the sample was deposited
    access_status - status of the sample ('public' or not)
    study - study title
    study_accession - project accession
    micro_environment - habitat
    macro_environment - high level habitat
    inserts_filtered - number of filtered inserts
    assembly_total_length - total assembled base pairs
    smORFs - total number of small ORFs
    ampsphere_amps - number of non-redundant AMPs
    is_host_associated - True or False 

