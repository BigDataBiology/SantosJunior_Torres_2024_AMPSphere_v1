**output_clustering_significance_levelIII.tsv**

**Description:**	Table from AMPSphere v2022-03 version. It consists in a table generated during
                        the validation of the clustering procedure. The validation consisted in the
                        alignment of 1000 randomly selected AMP candidates against the respective
                        representatives of their clusters.

The alignment was performed using Smith-Waterman algorithm in triplicate and
the identity/significance were calculated as in blast. Its columns are:
    
    query - randomly AMP candidate (excluded cluster representatives)
    target - cluster representative
    identity - identity calculated counting gaps
    gap_identity - identity calculated excluding gaps
    coverage - query coverage in percent of residues
    score - alignment score
    aln_len - alignment length
    evalue - expected-value from the exponential distribution by Altschull et al.
    sig. - '*' if e-value < 1e-5 else 'n.s.'
    family - AMP cluster accession at 3rd level (100-85-75% identity)
    replicate - replicate batch (1, 2, or 3)


**MD5 SUM:**	81d5ab8c604a785bbea3537532f5daef

**Size (MBytes):**	0.1034088134765625

**Content sample (first 5 items):**

query	target	identity	gap_identity	coverage	score	aln_len	evalue	sig.	family	replicate
AMP10.000_449	AMP10.074_887	96.7741935483871	96.7741935483871	100.0	173.0	31	5.464952928981811e-27	*	SPHERE-III.018_538	1
AMP10.000_799	AMP10.711_054	97.5	97.5	100.0	198.0	40	8.610889198974865e-32	*	SPHERE-III.029_362	1
AMP10.001_188	AMP10.123_547	46.666666666666664	48.275862068965516	100.0	70.0	30	8.936305862087231e-07	*	SPHERE-III.001_863	1
AMP10.001_319	AMP10.219_850	97.14285714285714	97.14285714285714	100.0	190.0	35	2.8137907983261325e-30	*	SPHERE-III.083_018	1
[...]
