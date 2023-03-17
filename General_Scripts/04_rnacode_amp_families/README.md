# RNAcode protein-coding sequences identification

This analysis used the software RNAcode to detect AMP genes with evolutionary 
signals of protein-coding sequences. To that, we used the clusters at the third 
level (100-85-75% id) containing at least 3 sequences, and assigned all its members
to the result of this analysis. Large clusters, with more than 500 members, were divided
into chunks of 100 randomly selected members and tested separed, the result obtained in 
2 out of the 3 chunks was considered the final verdict for the entire cluster -- most of
the cases it was consistent (all results were concordant).

The set of analysis include:

  (a) selection of clusters of at least 3 different AMP genes
  (b) subsampling of large clusters into 3 chunks of 100 AMPs
  (c) aligns gene sequences from the same cluster
  (d) RNAcode analysis
  (e) results compilation

To reproduce the analysis:

```
# to reproduce the RNAcode testing:
    $ python3 main.py
```

The RNAcode analysis outputs the following files to the *analysis/* folder:

*alignments/*

Folder containing clustal formatted alignments of genes from small clusters
with at least 3 sequences and less than 500. Alignments were made using MAFFT.

*rnacode_out/*

Results of RNAcode per alignment per chunk of each small cluster (SPHERE).
Format fixed by RNAcode software.

*large_clusters/aln/*

Folder containing clustal formatted alignments of genes from 3 chunks of each
large clusters (SPHERES) with at least 500 sequences.
Alignments were made using MAFFT.

*large_clusters/rnacode/*

Results of RNAcode per alignment per chunk of each large cluster (SPHERE).
Format fixed by RNAcode software.

*rnacode_output.tsv*

This table contains the significant protein-coding predicted segments
in the small clusters.
    
|**Column**|**Description**|
| :---: | :---: |
| HSS#	| Unique running number for each high scoring segment predicted in one RNAcode call |
| sense | Sense or anti-sense |
| frame | Translation Frame (1, 2, 3) |
| length | length of the predicted region in amino acids |
| start_aa | residue of the start of the predicted region |
| end_aa | residue of the end of the predicted region |
| rep_gene | representative gene (longest) |
| start | genomic coordinates (default 1) |
| end | end in genomic coordinates (arbitrary) |
| score | score of predicted region |
| p | p-value calculated |
| cluster | SPHERE code |
| amps | AMP sequences in the cluster |
| genes | Genes in the cluster |

*rnacode_out.merged.tsv*

This table contains the significant protein-coding predicted segments
in the small clusters and large clusters.
    
|**Column**|**Description**|
| :---: | :---: |
| HSS#	| Unique running number for each high scoring segment predicted in one RNAcode call |
| sense | Sense or anti-sense |
| frame | Translation Frame (1, 2, 3) |
| length | length of the predicted region in amino acids |
| start_aa | residue of the start of the predicted region |
| end_aa | residue of the end of the predicted region |
| rep_gene | representative gene (longest) |
| start | genomic coordinates (default 1) |
| end | end in genomic coordinates (arbitrary) |
| score | score of predicted region |
| p | p-value calculated |
| cluster | SPHERE code |
| amps | AMP sequences in the cluster |
| genes | Genes in the cluster |

*RNAcode_out_wlfam.tsv*

Two-column table containing the AMP access code and the decision if it is 
passing the RNAcode test (protein-coding sequence), Failling it (non-protein coding),
or could not be tested (not-tested).

The scripts in *utils/* folder are explained below:

	- genes_to_clusters.py:  **Function.** Select clusters of at least 3 genes and less than 500, align them and return them in clustal format 
	                         **Inputs.** AMPSphere_v.2022-03.fna.xz, SPHERE_v.2022-03.levels_assessment.tsv.gz
	                         **Outputs.** returns alignment files in clustal format for each cluster tested 

	- rnacode_call.py:       **Function.** Perform RNAcode test in the alignment files per cluster
	                         **Inputs.**   None
	                         **Outputs.**  generates two folders: analysis/alignments (containing the clustal format alignments per SPHERE),
	                                       and the analysis/rnacode_out/ (containing one file per SPHERE with the output from RNAcode)

	- largefams.py:          **Function.** Sample 3 times 100 sequences from the large SPHERES (with more than 500 genes) and perform the alignment
	                         **Inputs.** AMPSphere_v.2022-03.fna.xz, SPHERE_v.2022-03.levels_assessment.tsv.gz
	                         **Outputs.** returns alignment files in clustal format for each cluster tested 

	- rnacode_lfam.py:        **Function.** Perform RNAcode test in the alignment files per chunk of a large cluster and count the votes for True/False in each of them
	                         **Inputs.** None
	                         **Outputs.** generates two folders: analysis/large_clusters/aln/ (containing the clustal format alignments per SPHERE),
                                              and the analysis/large_clusters/rnacode/ (containing one file per SPHERE with the output from RNAcode)
                                              {ofile}')

	- lfam_arrange.py:       **Function.** Takes the results from large and small clusters and generate a final table with all results
	                         **Inputs.** rnacode outputs from both small and large clusters in the analysis/ folder
	                         **Outputs.** rnacode_out.merged.tsv, RNAcode_out_wlfam.tsv

Detailed inputs list:

| **Input file** | **Description** |
| :---: | :---: |
| AMPSphere_v.2022-03.fna.xz | AMPSphere nucleotide sequences |
| SPHERE_v.2022-03.levels_assessment.tsv.gz | Clustering key for the AMPSphere resource |

