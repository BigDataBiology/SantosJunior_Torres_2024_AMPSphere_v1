**families_all.count_core.tsv.gz**

**Description:**	TSV-file containing the data about the conservation of AMP families per species cluster in the
                        ProGenomes2 resource. To this analysis, the families were considered clusters of at least 8
                        AMP members. It was used to classify the AMP families into core (present in 95% or more 
                        of the genomes in a cluster), shell (present in less than 95% but in 50% or more of the
                        genomes in a cluster), or accessory (present in less than 50% of the genomes from a 
                        cluster). Its columns are as follows:

    families - AMP family accession (SPHERE level III) in the AMPSphere resource
    specI - species cluster in the ProGenomes2 resource
    counts - number of genomes from the given cluster containing the given AMP family
    total - total number of genomes composing the species cluster in question
    proportion - fraction of genomes in the cluster containing the given AMP family
    classification - classification of the given AMP family into core, shell or accessory

**MD5 SUM:**	97155d1118d04d16c5aea25e2179d926

**Size (MBytes):**	4.1961669921875e-05

**Content sample (first 5 items):**

families	specI	counts	total	proportion	classification
SPHERE-III.000_001	specI_v3_Cluster1726	1	26	3.8461538461538463	accessory
SPHERE-III.000_001	specI_v3_Cluster4808	3	15	20.0	accessory
SPHERE-III.000_001	specI_v3_Cluster95	8	4463	0.17925162446784673	accessory
SPHERE-III.000_002	specI_v3_Cluster2972	17	19	89.47368421052632	shell

