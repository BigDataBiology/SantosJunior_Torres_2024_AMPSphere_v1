**adjust_significant_function.csv**

**Description:**	TSV-file containing the data about the conservation of AMPs per species cluster in the
                        ProGenomes2 resource. It was used to classify the AMPs into core (present in 95% or more 
                        of the genomes in a cluster), shell (present in less than 95% but in 50% or more of the
                        genomes in a cluster), or accessory (present in less than 50% of the genomes from a 
                        cluster). Its columns are as follows:

    amp - AMP accession in the AMPSphere resource
    specI - species cluster in the ProGenomes2 resource
    counts - number of genomes from the given cluster containing the given AMP
    total - total number of genomes composing the species cluster in question
    proportion - fraction of genomes in the cluster containing the given AMP
    classification - classification of the given AMP into core, shell or accessory

**MD5 SUM:**	67ee29b19abb1f1479581c8be0c21bd8

**Size (MBytes):**	4.9591064453125e-05

**Content sample (first 5 items):**

amp	specI	counts	total	proportion	classification
AMP10.000_000	specI_v3_Cluster2367	12	26	46.15384615384615	accessory
AMP10.000_002	specI_v3_Cluster2367	13	26	50.0	shell
AMP10.000_004	specI_v3_Cluster95	3995	4463	89.51377996863097	shell
AMP10.000_005	specI_v3_Cluster855	14	16	87.5	shell
[...]
