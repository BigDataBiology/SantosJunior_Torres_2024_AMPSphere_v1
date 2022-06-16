# Core.py

Internal script
Run:

```
python3 uscripts/core.py
```

### : Description :

It calculates the proportion of core, shell and accessory AMPs and families.
Also including the high-quality sets curated from the assessment using the 5 tests 
from before.


### : Inputs :

1. **data/pgenomes_samples.tsv.gz:**

Resource from ProGenomes2 database, it consists in the
samples and their respective specI clusters, e.g.:

```
specI_v3_Cluster1	1049759.SAMN02436435
specI_v3_Cluster2	1001589.SAMN02436335
```

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

3. **data/SPHERE_v.2022-03.levels_assessment.tsv.gz:**

Table from AMPSphere v2022-03 version. It consists in a table with the clustering 
information in three different levels (100%, 100-85%, 100-85-75%) after conversion 
using an 8-letter alphabet. Its columns are:

    AMP accession - accession in AMPSphere v2022-03 resource
    evaluation vs. representative - summary of clustering procedure
    SPHERE_fam level I - cluster accession at 1st level (100% identity)
    SPHERE_fam level II - cluster accession at 2nd level (100-85% identity)
    SPHERE_fam level III - cluster accession at 3rd level (100-85-75% identity)

4. **data/quality_families.txt.gz:**

Table generated during the quality assessment. It contains the number of clusters
passing all the tests and which percent it represents for the cluster, also
contains the information about experimental evidence of at least 1 candidate in the
cluster. Its columns are:

    family - AMP cluster accession at 3rd level (100-85-75% identity)
    experimental_evidence - True or False(*)
    quality_candidates - number of peptides passing in all quality tests (**)
    total - total number of peptides in the cluster
    perc - percent of candidates passing in all tests

*  regarding to the identification of the AMP in metaproteomes or metatranscriptomes
   at least 1 candidate needed to become True

** simulatenously in RNAcode, AntiFam, and Terminal placement

5. **data/quality_candidates.txt.gz:**

List of candidates passing simultaneously in all quality tests (RNAcode, AntiFam, and
Terminal Placement), but without experimental evidence of translation or transcription.
This list is generated during the quality assessment and can be reobtained by using
the file ```quality_assessment.tsv.gz```.

6. **data/high_quality_candidates.txt.gz:**

List of candidates passing simultaneously in all quality tests (RNAcode, AntiFam, and
Terminal Placement) with experimental evidence of translation or transcription.
This list is generated during the quality assessment and can be reobtained by using
the file ```quality_assessment.tsv.gz```.


### : Outputs :

1. **families_all.count_core.tsv.gz:**

It is a table generated during the run of ```core.py``` script, and consists in the
analysis of families happening the samples from ProGenomes2 that were used in the
AMPSphere v2022-03. These families of AMPs were checked for each specI clusters in
AMPSphere and their prevalence in ProGenomes v2 from that given specI were accounted for.
Its columns are:

    family - AMP cluster accession at 3rd level (100-85-75% identity)
    specI - species cluster from ProGenomes2
    counts - number of genomes from a given specI in which we found the family
    total - number of genomes from a given specI
    proportion - percent of genomes from a specI containing a given family
    classification - classification into core, shell and accessory (***)
    
*** - core: >90% of prevalence,
      shell: >50% and <90% of prevalence,
      accessory: < 50% of prevalence in the genomes from a given specI

2. **amps_all.count_core.tsv.gz:**

It is a table generated during the run of ```core.py``` script, and consists in the
analysis of AMP candidates happening the samples from ProGenomes2 that were used in the
AMPSphere v2022-03. These AMPs were checked for each specI clusters in AMPSphere and their
prevalence in ProGenomes v2 from that given specI were accounted for. Its columns are:

    amp - accession of AMP in AMPSphere v2022-03
    specI - species cluster from ProGenomes2
    counts - number of genomes from a given specI in which we found the AMP
    total - number of genomes from a given specI
    proportion - percent of genomes from a specI containing a given AMP
    classification - classification into core, shell and accessory (***)
    
*** - core: >90% of prevalence,
      shell: >50% and <90% of prevalence,
      accessory: < 50% of prevalence in the genomes from a given specI

3. **core_analysis.svg:**

Stacked Bar chart oriented horizontally containing the three classes: core, shell
and accessory AMPs and families, high-quality or not.

