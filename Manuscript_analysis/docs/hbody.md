# Hbody.py

Internal script
Run:

```
python3 uscripts/hbody.py
```

### : Description :

It generates a heatmap with the overlap of AMP candidates from 
differen human body sites. It calculates the overlap in counts and
in % by the column maximum.

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
    
### : Outputs :

1. **percent_bycolumn_overlap_hbody.tsv**

Squared table containing human body sites. The values are the of column AMPs
shared with the row. 

2. **figure_1d_alt_heatmap_humanbodysites_overlap.svg**

Heatmap contaning the different human body sites with the color scale
representing the percent of AMPs shared from the row with the column.
E.g.: 19% of AMPs from mouth (column) are also present in the gut (row).

