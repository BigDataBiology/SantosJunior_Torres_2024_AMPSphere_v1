## Genome context

We used a similar approach to that described in Río et al. (2022).
The data in here was obtained running the pipeline from GeccoViz.

The outputs are in the `data/` folder, and consist of:

| **Output file** | **Description** |
| :---: | :---: |
| AMPs_hit_genomes.txt.xz | output list of AMPs with hits in the curated set of genomes |
| Genomic_context_frequency_KEGG_pathways.tsv.xz | output table bringing the genome context frequency for each KEGG pathway |
| hits.folded.tab.xz | output table containing the AMP hit number and genomes |
| neigh_VCS_higher_0.9_3_or_more_hits.tab.xz | output table containing the hits for AMPs and their gene neighbors conserved in the genome context |
| permutation_tests.tab.xz | output table containing the 10,000 permutations of sampling 50,000 full-length protein families of different sizes and counting the number of them with conserved genome contexts and genome contexts containing specific genes/terms |
| strong_hits.txt.xz | output list of AMPs with hits characterizing them as strong homologs (high identity and significance) in the curated set of genomes |

____

### References

Álvaro Rodríguez del Río, Joaquín Giner-Lamia, Carlos P. Cantalapiedra, Jorge Botas, Ziqi Deng,
Ana Hernández-Plaza, Lucas Paoli, Thomas S.B. Schmidt, Shinichi Sunagawa, Peer Bork, Luis Pedro Coelho,
Jaime Huerta-Cepas. Functional and evolutionary significance of unknown genes from uncultivated taxa.
bioRxiv 2022.01.26.477801; doi: https://doi.org/10.1101/2022.01.26.477801 
