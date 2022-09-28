**dramp_v3.features.tsv**

**Description:**	The AMPs from [DRAMP v3](DRAMP.cpu-bioinfor.org/) were filtered to eliminate the peptides
                        containing non-cannonical residues. 

                        This table contains features calculated with [Macrel](macrel.readthedocs.io/) pipeline.
                        Its columns are:
    
    Access - AMP accession in DRAMP v3
    sequence - amino acids sequence
    group - decoy column, label 'Unk'
    tinyAA - proportion of (A+C+G+S+T) in the peptide
    smallAA - proportion of (A+B+C+D+G+N+P+S+T+V) in the peptide
    aliphaticAA - proportion of (A+I+L+V) in the peptide
    aromaticAA - proportion of (F+H+W+Y) in the peptide
    nonpolarAA - proportion of (A+C+F+G+I+L+M+P+V+W+Y) in the peptide
    polarAA - proportion of (D+E+H+K+N+Q+R+S+T+Z) in the peptide
    chargedAA - proportion of (B+D+E+H+K+R+Z) in the peptide
    basicAA - proportion of (H+K+R) in the peptide
    acidicAA - proportion of (B+D+E+Z) in the peptide
    charge - calculated charge at pH 7 with pKscale "EMBOSS"
    pI - isoelectric point using pKscale "EMBOSS"
    aindex - relative volume occupied by aliphatic side chains (AVIL)
    instaindex - stability of a protein based on its amino acid composition
    boman - overall estimate of the potential of a peptide to bind to membranes or other
            proteins as receptor
    hydrophobicity - average hydrophobicity using scale "KyteDoolittle"
    hmoment - quantitative measure of the amphiphilicity perpendicular to the axis of any
              periodic peptide structure, such as the alpha-helix or beta-sheet
    SA.Group1.residue0 - CTDD at first residue, using solvent accessibility, class I
    SA.Group2.residue0 - CTDD at first residue, using solvent accessibility, class II
    SA.Group3.residue0 - CTDD at first residue, using solvent accessibility, class III
    HB.Group1.residue0 - CTDD at first residue, using FET from water to lipid phase, CI
    HB.Group2.residue0 - CTDD at first residue, using FET from water to lipid phase, CII
    HB.Group3.residue0 - CTDD at first residue, using FET from water to lipid phase, CIII

Solvent Accessibility (SA) was adopted as in previous studies (Dubchak et al. 1995, 1999),
however, the new feature HB was adapted from Von Heijne and Blomberg, 1979. Classes adopted
to the sequence encoding of the distribution at the first residue of each class:

 | **Properties** | **Class I** | **Class II** | **Class III** |
 | :---: | :---: | :---: | :---: |
 | SA | A, L, F, C, G, I, V, W | R, K, Q, E, N, D | M, S, P, T, H, Y |
 | FT | I,L,V,W,A,M,G,T | F,Y,S,Q,C,N | P,H,K,E,D,R |


**MD5 SUM:**	db93deefeae4304941ab3891d822eaaf

**Size (MBytes):**	1.375802993774414

**Content sample (first 5 items):**

Access	sequence	group	tinyAA	smallAA	aliphaticAA	aromaticAA	nonpolarAA	polarAA	chargedAA	basicAA	acidicAA	charge	pI	aindex	instaindex	boman	hydrophobicity	hmoment	SA.Group1.residue0	SA.Group2.residue0	SA.Group3.residue0	HB.Group1.residue0	HB.Group2.residue0	HB.Group3.residue0
DRAMP00001_GENERAL	GSGVIPTISHECHMNSFQFVFTCCS	Unk	0.44	0.6	0.16	0.2	0.56	0.44	0.12	0.08	0.04	-0.6343022196857381	6.47406005859375	54.4	36.047999999999995	0.4847999999999999	0.2836	0.29962104127016576	4.0	44.0	8.0	4.0	8.0	24.0
DRAMP00002_GENERAL	WKSESVCTPGCVTGLLQTCFLQTITCNCKISK	Unk	0.46875	0.59375	0.21875	0.0625	0.53125	0.46875	0.125	0.09375	0.03125	1.8234093809336598	8.21258544921875	79.0625	41.071875000000006	0.5612499999999999	0.12874999999999998	0.3779634686290898	3.125	6.25	9.375	3.125	9.375	6.25
DRAMP00003_GENERAL	ADRGWIKTLTKDCPNVISSICAGTIITACKNCA	Unk	0.48484848484848486	0.6666666666666666	0.3333333333333333	0.030303030303030304	0.5757575757575758	0.42424242424242425	0.18181818181818182	0.12121212121212122	0.06060606060606061	1.8543897017426985	8.33978271484375	91.8181818181818	24.975757575757576	0.886060606060606	0.10969696969696967	0.565095576997773	3.0303030303030303	6.0606060606060606	24.242424242424242	3.0303030303030303	39.39393939393939	6.0606060606060606
DRAMP00004_GENERAL	NNTIKDFDLDLKTNKKDTATPYVGSRYLCTPGSCWKLVCFTTTVK	Unk	0.35555555555555557	0.6222222222222222	0.2	0.1111111111111111	0.4666666666666667	0.5333333333333333	0.24444444444444444	0.15555555555555556	0.08888888888888889	2.8845677391613274	8.89764404296875	64.88888888888889	21.70888888888889	1.7359999999999998	-0.06599999999999992	0.46406919588308515	8.88888888888889	2.2222222222222223	6.666666666666667	6.666666666666667	2.2222222222222223	11.11111111111111
[...]
