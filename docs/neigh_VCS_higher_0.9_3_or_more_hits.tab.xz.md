**neigh_VCS_higher_0.9_3_or_more_hits.tab**

**Description:**	TSV-file containing families (after extending with our genomes) with 3 or more members
                        with conserved context (ie. conservation score > 0.9 with known function genes). These
                        are good candidates for identifying novel families potentially involved in known processes.

| **columns** | **description** |
| :---: | :---: |
| AMP | AMP name |
| db | db to which the functional term in "funct_name" belongs (eggnog, kegg, etc) |
| funct_name | name of the gene family with which the AMP has conserved genomic context. |
| position | position of the known function gene relative to the novel family (ie. +1 means that the neighbour gene is the next gene to the novel family) |
| conservation score | number of members of the novel family which, in the position indicated, have a neighbour with the function indicated, divided by the number of genes the novel family has (the closer to 1 the better) |
| % genes in different strand in that position | % genes, in the positioin indicated and with the function indicated, that are in the contrary strand to the novel family members (the closer to 0 the better) |
| % genes in different strand in between the target gene and the conserved gene | in the case the neighbour gene is not in the +-1 position, proportion of genes which are in the contrary strand in between the novel family member and the gene un the position indicated (the closer to 0 the better ) |
| % high distance genes in between the target gene and the conserved function gene | proportion of genes in between the novel family members and the gene in the position indicated & with the function indicated which are separated more than 100nt away in the direction of the novel family members (the closer to 0 the better ) |
| description | description of gene/pathway function |


**MD5 SUM:**	7bf6d5d3eb45a275df097544c5fd382b

**Size (MBytes):**	1.789968

**Content sample (first 5 items):**

AMP	db	funct_name	position	v_score	% genes in different strand in that position	% genes in different strand in between the target gene and the conserved gene	% high distance genes in between the target gene and the conserved function gene	description
AMP10.000_003	ogs	COG0806	1	0.9956896551724138	0.0	0.0	0.7316017316017316	ribosome binding
AMP10.000_003	ogs	1V6HD	1	0.9956896551724138	0.0	0.0	0.7316017316017316	
AMP10.000_003	ogs	24I1G	1	0.9956896551724138	0.0	0.0	0.7316017316017316	
AMP10.000_003	ogs	3WIMH	1	0.9956896551724138	0.0	0.0	0.7316017316017316
[...]
