**output_clustering_significance_levelII.tsv**

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
    family - AMP cluster accession at 2nd level (100-85% identity)
    replicate - replicate batch (1, 2, or 3)


**MD5 SUM:**	b8ec7ec448e2525c6b355651e006adf2

**Size (MBytes):**	0.09694671630859375

**Content sample (first 5 items):**

query	target	identity	gap_identity	coverage	score	aln_len	evalue	sig.	family	replicate
AMP10.000_069	AMP10.308_984	63.63636363636363	100.0	100.0	103.5	33	1.6874935210442608e-13	*	SPHERE-II.089_007	1
AMP10.000_108	AMP10.473_048	60.714285714285715	70.83333333333333	100.0	76.5	28	3.904178998690963e-08	*	SPHERE-II.001_068	1
AMP10.000_198	AMP10.015_876	96.7741935483871	96.7741935483871	100.0	165.0	31	2.0409021195010538e-25	*	SPHERE-II.024_205	1
AMP10.000_456	AMP10.017_849	96.66666666666667	96.66666666666667	100.0	150.0	30	1.7519583240735643e-22	*	SPHERE-II.019_933	1
[...]
