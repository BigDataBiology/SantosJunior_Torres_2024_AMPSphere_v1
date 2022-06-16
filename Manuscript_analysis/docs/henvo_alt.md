# Henvo_alt.py

Internal script
Run:

```
python3 uscripts/henvo_alt.py
```

### : Description :

It generates a heatmap with the overlap of AMP candidates from 
different high level environments.
It calculates the overlap in counts and in % by the column maximum.

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

1. **percent_bycolumn_overlap_higherenvo.tsv**

Squared table containing high level environments.
The values are number of AMPs from the environment in
the row are shared with the environment in the column. 

2. **figure_1c_alt_heatmap_environments_overlap.svg**

Heatmap contaning the different high level environments with the color scale
representing the percent of AMPs shared from the row with the column.
E.g.: 19% of AMPs from human guts (column) are also present in the plant/soil (row).

3. **lowlevel_habitat.tsv**

Squared table containing all environments with at least 100 AMP candidates.
The values are number of AMPs from the environment in
the row are shared with the environment in the column. 

4. **low_level_habitats_ge100pep.svg**

Heatmap contaning the different environments with at least 100 AMPs, the color scale
represent the percent of AMPs shared from the row with the column.
E.g.: 19% of AMPs from human guts (column) are also present in the plant-associated (row).

