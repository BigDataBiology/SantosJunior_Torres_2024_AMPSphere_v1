# Sup_table_high_environments.py

Internal scriptquality_families.txt.gz
Run:

```
python3 uscripts/sup_table_high_environments.py
```

### : Description :

It generates the supplementary table with data for the higher environments. 

### : Inputs :

1. **data/SPHERE_v.2022-03.levels_assessment.tsv.gz:**

Table from AMPSphere v2022-03 version. It consists in a table with the clustering 
information in three different levels (100%, 100-85%, 100-85-75%) after conversion 
using an 8-letter alphabet. Its columns are:

    AMP accession - accession in AMPSphere v2022-03 resource
    evaluation vs. representative - summary of clustering procedure
    SPHERE_fam level I - cluster accession at 1st level (100% identity)
    SPHERE_fam level II - cluster accession at 2nd level (100-85% identity)
    SPHERE_fam level III - cluster accession at 3rd level (100-85-75% identity)

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
    
### : Outputs :

1. **supplementary_table_S2.tsv**

Table bringing the following columns:

    high level environment - classification of habitat into high level
    habitats - habitat from which the samples originated
    samples - number of samples associated to a high level habitat
    redundant AMPs - number of redundant AMP genes in an environment
    non-redundant AMPs - number of non-redundant AMP genes in an environment
    AMP clusters - number of AMP clusters (include singletons)
    AMP families - number of AMP clusters with at least 8 peptides (families)


