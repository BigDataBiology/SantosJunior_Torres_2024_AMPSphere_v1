**output_clustering_significance_levelI.tsv**

**Description:**	Table from AMPSphere v2022-03 version. It consists in a table generated during
                        the validation of the clustering procedure. The validation consisted in the
                        alignment of 1000 randomly selected AMP candidates against the respective
                        representatives of their clusters in a given level.

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
    family - AMP cluster accession at 1st level (100% identity)
    replicate - replicate batch (1, 2, or 3)

**MD5 SUM:**	d76bb306c64c8bae181a527aa9f31b20

**Size (MBytes):**	0.08940887451171875

**Content sample (first 5 items):**

query	target	identity	gap_identity	coverage	score	aln_len	evalue	sig.	family	replicate
AMP10.000_247	AMP10.044_080	94.44444444444444	100.0	100.0	167.5	36	7.22130565500068e-26	*	SPHERE-I.009_727	1
AMP10.000_290	AMP10.249_586	97.67441860465117	100.0	100.0	211.0	43	2.519731147495874e-34	*	SPHERE-I.003_554	1
AMP10.000_448	AMP10.384_547	97.67441860465117	97.67441860465117	100.0	215.0	43	4.221390166631519e-35	*	SPHERE-I.004_353	1
AMP10.001_278	AMP10.671_642	97.43589743589743	97.43589743589743	100.0	193.0	39	8.066775323213292e-31	*	SPHERE-I.004_163	1
[...]
