# Graph_enrichment.py

Internal script
Run:

```
python3 uscripts/graph_enrichment.py
```

### : Description :

It plots the results of the mapping of AMP candidates against large proteins from GMGC
returning horizontal bar plots for the ortholog groups enrichment and proportion of 
COG classes among the homologs.


### : Inputs :

1. **data/adjust_significant_function.csv.xz**

Enrichment analysis of functions among those ortholog groups identified by mapping
AMPs against GMGC. This analysis was done using all proteins in GMGC not only the
unigenes, which were deduplicated prior counting. Its columns are as follow:

    eggnog_OG - EggNOG ortholog group
    count_AMP - number of AMPs matching to an OG
    count_GMGC - number of proteins from GMGC matching an OG
    total_AMP - total AMPs matched against GMGC
    total_GMGC - total proteins in GMGC
    p_value - p_value calculated from the hypergeometric space
    amp_fraction - proportion of total AMPs matching an OG
    GMGC_fraction - proportion of large proteins matching an OG
    times - fold enrichment (AMP fraction / GMGC fraction)
    p_adjust - p-value adjusted using Bonferroni
    sig - if p<0.05
    label - decoy

2. **data/amp_COG.tsv.xz**

Table contaning the COG classes for the mapped AMPs. It does contain about 6k entries
because most of the OGs from GMGC were not in the metadata table relative to the 
EggNOG mapping. Its columns are:

    cog_class - class from [COG database](https://ecoliwiki.org/colipedia/index.php/Clusters_of_Orthologous_Groups_(COGs))
    count_AMP - mapped AMPs
    fraction - proportion of mapped AMPs


### : Outputs :

1. **top_og_enrichment.svg**

Horizontal bar plot with the `x axis` as the Enrichment (fold) and the `y axis` as the EggNOG ortholog groups.

2. **most_frequent_og.svg**

Horizontal bar plot with the `x axis` as the AMP orthologs and the `y axis` as the EggNOG ortholog groups.
    
3. **most_frequent_cogclass.svg**

Horizontal bar plot with the `x axis` as the AMP orthologs and the `y axis` as the COG classes.


