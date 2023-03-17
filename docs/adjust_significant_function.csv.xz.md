**adjust_significant_function.csv**

**Description:**	Enrichment analysis of functions among those ortholog groups identified by mapping
                        AMPs against GMGC. This analysis was done using all proteins in GMGC not only the
                        unigenes, which were deduplicated prior counting. Its columns are as follows:

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

**MD5 SUM:**	74e01af97c6ededa2a098cf5c7dda4f6

**Size (MBytes):**	0.7440414428710938

**Content sample (first 5 items):**

"eggnog_OG"	"count_AMP"	"count_GMGC"	"total_AMP"	"total_GMGC"	"p_value"	"amp_fraction"	"GMGC_fraction"	"times"	"p_adjust"
"COG3547"	1356	3326230	156711	9180087363	0	0.0086528705706683	0.000362330974474845	23.8811230069678	0
"1VKBK"	244	4072	156711	9180087363	0	0.00155700620888132	4.4356876345339e-07	3510.18001529074	0
"24UUY"	243	4064	156711	9180087363	0	0.0015506250358941	4.42697312051713e-07	3502.67551593574	0
"2C6KI"	244	4072	156711	9180087363	0	0.00155700620888132	4.4356876345339e-07	3510.18001529074	0
[...]
