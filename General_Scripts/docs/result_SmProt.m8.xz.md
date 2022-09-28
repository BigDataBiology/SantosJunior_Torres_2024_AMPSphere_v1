**result_SmProt.m8**

**Description:**	MMSeqs2 easy-search results using candidate AMPs from AMPSphere against
                        the [SmProt 2 database](http://bigdata.ibp.ac.cn/SmProt/) after quality
                        filtering leaving only mhits with identities above 75% covering at least
                        50% of the shorter sequence and with e-value < 1e-5.

| **column position** | **description** |
| :---: | :---: |
| 1 | AMP access code from the AMPSphere resource |
| 2 | AMP access code from the SmProt 2 database |
| 3 | alignment E-value |
| 4 | number of gaps openned in the alignment |
| 5 | Percentage of identical matches |
| 6 | Number of identical matches |
| 7 | position of the alignment start in the query sequence |
| 8 | position of the alignment stop in the query sequence  |
| 9 | query total length in residues |
| 10 | position of the alignment start in the target sequence |
| 11 | position of the alignment stop in the target sequence |
| 12 | target total length in residues |
| 13 | alignment length in residues |
| 14 | raw alignment score |
| 15 | BitScore calculated from the alignment |
| 16 | Alignment as string. Each position contains either M (match), D (deletion, gap in query), or I (Insertion, gap in target) |
| 17 | query sequence |
| 18 | target sequence |
| 19 | Header of Query sequence |
| 20 | Header of Target sequence |
| 21 | Aligned query sequence with gaps |
| 22 | Aligned target sequence with gaps |
| 23 | frame of alignment with query, null for proteins |
| 24 | frame of alignment with target, null for proteins |
| 25 | number of mismatches in the alignment |
| 26 | query coverage % |
| 27 | target coverage % |


**MD5 SUM:**	3a43421db3bf220c9b08a3dcd092c648

**Size (MBytes):**	4.435081481933594

**Content sample (first 5 items):**

AMP10.073_889	SPROECO6162	7.714E-07	0	78.500	22	1	28	28	39	66	68	28	100	44	28M	LFKNGAVEATKVGALSKSQLAAFLDANV	LDEIADEYQGKLTVAKLNIDQNPGTAPKYGIRGIPTLLLFKNGEVAATKVGALSKGQLKEFLDANLA*	AMP10.073_889 | SPHERE-III.000_437	SPROECO6162	LFKNGAVEATKVGALSKSQLAAFLDANV	LFKNGEVAATKVGALSKGQLKEFLDANL			6	1.000	0.412
AMP10.073_889	SPROECO6162	7.714E-07	0	78.500	22	1	28	28	39	66	68	28	100	44	28M	LFKNGAVEATKVGALSKSQLAAFLDANV	LDEIADEYQGKLTVAKLNIDQNPGTAPKYGIRGIPTLLLFKNGEVAATKVGALSKGQLKEFLDANLA*	AMP10.073_889 | SPHERE-III.000_437	SPROECO6162	LFKNGAVEATKVGALSKSQLAAFLDANV	LFKNGEVAATKVGALSKGQLKEFLDANL			6	1.000	0.412
AMP10.073_889	SPROECO281	7.714E-07	0	78.500	22	1	28	28	43	70	72	28	100	44	28M	LFKNGAVEATKVGALSKSQLAAFLDANV	IAPILDEIADEYQGKLTVAKLNIDQNPGTAPKYGIRGIPTLLLFKNGEVAATKVGALSKGQLKEFLDANLA*	AMP10.073_889 | SPHERE-III.000_437	SPROECO281	LFKNGAVEATKVGALSKSQLAAFLDANV	LFKNGEVAATKVGALSKGQLKEFLDANL			6	1.000	0.389
AMP10.073_889	SPROECO281	7.714E-07	0	78.500	22	1	28	28	43	70	72	28	100	44	28M	LFKNGAVEATKVGALSKSQLAAFLDANV	IAPILDEIADEYQGKLTVAKLNIDQNPGTAPKYGIRGIPTLLLFKNGEVAATKVGALSKGQLKEFLDANLA*	AMP10.073_889 | SPHERE-III.000_437	SPROECO281	LFKNGAVEATKVGALSKSQLAAFLDANV	LFKNGEVAATKVGALSKGQLKEFLDANL			6	1.000	0.389
AMP10.073_889	SPROECO10269	7.714E-07	0	78.500	22	1	28	28	44	71	73	28	100	44	28M	LFKNGAVEATKVGALSKSQLAAFLDANV	MIAPILDEIADEYQGKLTVAKLNIDQNPGTAPKYGIRGIPTLLLFKNGEVAATKVGALSKGQLKEFLDANLA*	AMP10.073_889 | SPHERE-III.000_437	SPROECO10269	LFKNGAVEATKVGALSKSQLAAFLDANV	LFKNGEVAATKVGALSKGQLKEFLDANL			6	1.000	0.384
[...]
