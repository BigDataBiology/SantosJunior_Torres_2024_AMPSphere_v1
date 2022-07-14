# Family_sizes.py

Internal script
Run:

```
python3 uscripts/family_sizes.py
```

### : Description :

It generates two histograms with the distribution of families by the number of 
genera they affiliate. One graph contains the general view and another only the
families with 20 genera or more.

### : Inputs :

1. **data/SPHERE_v.2022-03.levels_assessment.tsv.gz**

Table from AMPSphere v2022-03 version. It consists in a table with the clustering 
information in three different levels (100%, 100-85%, 100-85-75%) after conversion 
using an 8-letter alphabet. Its columns are:

    AMP accession - accession in AMPSphere v2022-03 resource
    evaluation vs. representative - summary of clustering procedure
    SPHERE_fam level I - cluster accession at 1st level (100% identity)
    SPHERE_fam level II - cluster accession at 2nd level (100-85% identity)
    SPHERE_fam level III - cluster accession at 3rd level (100-85-75% identity)

2. **data/quality_families.txt.gz**

Table generated during the quality assessment. It contains the number of clusters
passing all the tests and which percent it represents for the cluster, also
contains the information about experimental evidence of at least 1 candidate in the
cluster. Its columns are:

    family - AMP cluster accession at 3rd level (100-85-75% identity)
    experimental_evidence - True or False('*')
    quality_candidates - number of peptides passing in all quality tests ('**')
    total - total number of peptides in the cluster
    perc - percent of candidates passing in all tests

'*'  regarding to the identification of the AMP in metaproteomes or metatranscriptomes
   at least 1 candidate needed to become True

'**' simulatenously in RNAcode, AntiFam, and Terminal placement

3. **data/taxonomy_annotation.tsv.gz**

Table from AMPSphere v2022-03 version. It consists in a table generated by filtering
the taxonomy assigned to contigs and AMP candidates. The AMP candidates assigned to
genus or species levels were considerd. This table has the info related to the
lowest taxonomic level each AMP had in AMPSphere. Its columns are:
    
    amp - accession in AMPSphere v2022-03 resource
    taxid - accession in NCBI taxID database
    level - lowest taxonomy classification (Unknown, Root, Kingdom, Phylum,
            Class, Order, Family, Genus, Species)
    source - taxon assigned to the contig from which the AMP was predicted
    fixed - one term designating taxonomy, here stands for genus
    
### : Outputs :

1. **qc_families_numberofgenera.svg**

Histogram of the number of families happening by different number of genera.
It shows in the `x axis` the number of genera in a given family and in the `y axis`
the number of families.

2. **qc_families_numberofgenera_ge100.svg**

Histogram of the number of families happening by different number of genera.
Only families contaning 100 or more genera made to this plot.
It shows in the `x axis` the number of genera in a given family and in the `y axis`
the number of families.
