**macrel_trainpos.features.tsv**

**Description:**	The AMPs from [Macrel traning set](https://github.com/BigDataBiology/macrel/tree/master/train)
                        were used to calculate features with [Macrel](macrel.readthedocs.io/) pipeline.
                        Its columns are:
    
    Access - AMP accession in Macrel training set
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

**MD5 SUM:**	f052ed94f0ce4b3d369ba17b72e7c448

**Size (MBytes):**	0.3280792236328125

**Content sample (first 5 items):**

Access	sequence	group	tinyAA	smallAA	aliphaticAA	aromaticAA	nonpolarAA	polarAA	chargedAA	basicAA	acidicAA	charge	pI	aindex	instaindex	boman	hydrophobicity	hmoment	SA.Group1.residue0	SA.Group2.residue0	SA.Group3.residue0	HB.Group1.residue0	HB.Group2.residue0	HB.Group3.residue0
AMP_0	AACSDRAHGHICESFKSFCKDSGRNGVKLRANCKKTCGLC	Unk	0.475	0.6	0.2	0.1	0.5	0.5	0.325	0.25	0.075	5.274523077950661	8.89642333984375	46.50000000000001	31.572500000000005	2.302249999999999	-0.1915	0.6259965040256069	2.5	12.5	10.0	2.5	7.5	12.5
AMP_1	AAEFPDFYDSEEQMGPHQEAEDEKDRADQRVLTEEEKKELENLAAMDLELQKIAEKFSQR	Unk	0.18333333333333332	0.35	0.23333333333333334	0.08333333333333333	0.38333333333333336	0.6166666666666667	0.4666666666666667	0.15	0.31666666666666665	-10.764340211740313	4.01617431640625	55.49999999999999	60.02683333333333	3.351	-0.30033333333333334	0.6514605750866599	1.6666666666666667	5.0	8.333333333333332	1.6666666666666667	6.666666666666667	5.0
AMP_2	AAFFAQQKGLPTQQQNQVSPKAVSMIVNLEGCVRNPYKCPADVWTNGVGNTHNVDKTKILTIDEVATDLRRNIKEAENCINTYFNGEKMNQGQYDAMVSLAFNVGCGNIKTYYSKTQGKRVATTIYRAAQAENWILMCNRIEDFNKSGGRVLKGLQNRRAKEKALCLGE	Unk	0.2958579881656805	0.5266272189349113	0.27218934911242604	0.08284023668639054	0.5088757396449705	0.4911242603550296	0.23076923076923078	0.14201183431952663	0.08875739644970414	8.041296197687904	9.45526123046875	73.31360946745562	33.66568047337277	1.9657396449704136	-0.08792899408284022	0.7730155738254879	0.591715976331361	3.5502958579881656	6.508875739644971	0.591715976331361	1.7751479289940828	4.733727810650888
AMP_3	AAFRGCWTKNYSPKPCL	Unk	0.4117647058823529	0.5882352941176471	0.17647058823529413	0.17647058823529413	0.6470588235294118	0.35294117647058826	0.17647058823529413	0.17647058823529413	0.0	2.9134739326759838	9.49700927734375	34.705882352941174	37.811764705882354	1.26	-0.039999999999999994	0.2215095489054466	5.88235294117647	23.52941176470588	47.05882352941176	5.88235294117647	17.647058823529413	23.52941176470588
[...]
