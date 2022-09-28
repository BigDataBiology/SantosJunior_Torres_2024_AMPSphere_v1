**ampsphere_v2022-03.features.tsv**

**Description:**	Features of the candidate AMPs in the AMPSphere resource calculated with the Macrel scripts.
                        It consists in a TSV-file with the following columns:
                        
| **columns** | **description** |
| :---: | :---: |
| Access | AMPSphere's c_AMP access code |
| sequence | peptide sequence |
| group | dummy variable introduced by macrel, the 'Unk' stands for unknown state prior prediction, which here does not mean anything |
| tinyAA | Percent of residues of tiny amino acids (A + C + G + S + T) |
| smallAA | Percent of residues of small amino acids (A + B + C + D + G + N + P + S + T + V) |
| aliphaticAA | Percent of residues of aliphatic amino acids (A + I + L + V) |
| aromaticAA | Percent of residues of aromatic amino acids (F + H + W + Y) |
| nonpolarAA | Percent of residues of non-polar amino acids (A + C + F + G + I + L + M + P + V + W + Y) |
| polarAA | Percent of residues of polar amino acids (D + E + H + K + N + Q + R + S + T + Z) |
| chargedAA | Percent of residues of charged amino acids (B + D + E + H + K + R + Z) |
| basicAA | Percent of residues of basic amino acids (H + K + R) |
| acidicAA | Percent of residues of acidic amino acids (B + D + E + Z) |
| charge | Peptide charge at pH 7.0 using "EMBOSS" pk-scale |
| pI | Peptide isoelectric point using "EMBOSS" pk-scale |
| aindex | Relative volume occupied by aliphatic side chains (A, V, I, and L) in the peptide chain |
| instaindex | Stability of a protein based on its dipeptide composition, it is used to determine whether the peptide will be stable in a test tube. |
| boman | Sum of the solubility values for all residues in a sequence, it might give an overall estimate of the potential of a peptide to bind to membranes or other proteins as receptors, to normalize it is divided by the number of residues |
| hydrophobicity | Peptide's hydrophobicity using "KyteDoolittle" scale |
| hmoment | Quantitative measure of the amphiphilicity perpendicular to the axis of any periodic peptide structure, such as the alpha-helix or beta-sheet using angle of 100º and a window of 11 residues |
| SA.Group1.residue0 | Distribution parameter (CTDD) at first occurrence in the sequence using amino acids  group 1 (A,L,F,C,G,I,V,W) clustered by solvent accessibility |
| SA.Group2.residue0 | Distribution parameter (CTDD) at first occurrence in the sequence using amino acids  group 2 (R,K,Q,E,D,N)  clustered by solvent accessibility |
| SA.Group3.residue0 | Distribution parameter (CTDD) at first occurrence in the sequence using amino acids  group 3 (M,S,P,T,H,Y)  clustered by solvent accessibility |
| HB.Group1.residue0 | Distribution parameter (CTDD) at first occurrence in the sequence using amino acids  group 1 (I,L,V,W,A,G,T,M)  clustered by the free energy to transfer from water to lipophilic phase |
| HB.Group2.residue0 | Distribution parameter (CTDD) at first occurrence in the sequence using amino acids  group 2 (F,Y,S,Q,C,N) clustered by the free energy to transfer from water to lipophilic phase |
| HB.Group3.residue0 | Distribution parameter (CTDD) at first occurrence in the sequence using amino acids  group 3 (P,H,K,E,D,R) clustered by the free energy to transfer from water to lipophilic phase |

**MD5 SUM:**	9fa290deb324579109bb3b31594daa56

**Size (MBytes):**	75.22193050384521

**Content sample (first 5 items):**

Access	sequence	group	tinyAA	smallAA	aliphaticAA	aromaticAA	nonpolarAA	polarAA	chargedAA	basicAA	acidicAA	charge	pI	aindex	instaindex	boman	hydrophobicity	hmoment	SA.Group1.residue0	SA.Group2.residue0	SA.Group3.residue0	HB.Group1.residue0	HB.Group2.residue0	HB.Group3.residue0
AMP10.000_000	KKVKSIFKKALAMMGENEVKAWGIGIK	Unk	0.25925925925925924	0.37037037037037035	0.3333333333333333	0.07407407407407407	0.5925925925925926	0.4074074074074074	0.3333333333333333	0.25925925925925924	0.07407407407407407	4.977300016786144	10.92608642578125	90.37037037037037	-18.348148148148148	0.6107407407407408	0.03703703703703704	0.6774924662026925	11.11111111111111	3.7037037037037033	18.51851851851852	11.11111111111111	18.51851851851852	3.7037037037037033
AMP10.000_001	FFGIGQQEMTLEEIGDKFGLTRERVRQIKEKAIRRLRQSNRSKLLKSYLG	Unk	0.22	0.28	0.24	0.08	0.44	0.56	0.36	0.24	0.12	5.981366836855313	11.25164794921875	85.8	59.397999999999996	2.963999999999999	-0.2836	0.7922628140721152	2.0	12.0	18.0	6.0	2.0	16.0
AMP10.000_002	KRVKSFFKGYMRAIEINAALMYGYRPK	Unk	0.2222222222222222	0.3333333333333333	0.25925925925925924	0.18518518518518517	0.6296296296296297	0.37037037037037035	0.2962962962962963	0.25925925925925924	0.037037037037037035	5.974127486750387	10.93487548828125	65.18518518518519	49.02222222222221	1.7577777777777774	-0.11148148148148142	0.6922105739402304	11.11111111111111	3.7037037037037033	18.51851851851852	11.11111111111111	18.51851851851852	3.7037037037037033
AMP10.000_003	GRVIGKQGRIAKAIRVVMRAAAVRVDEKVLVEID	Unk	0.23529411764705882	0.5	0.5	0.0	0.6176470588235294	0.38235294117647056	0.35294117647058826	0.23529411764705882	0.11764705882352941	3.9795054578216305	11.52386474609375	131.76470588235293	10.691176470588236	1.786764705882353	-0.057647058823529385	0.8805951799204162	2.941176470588235	5.88235294117647	52.94117647058824	2.941176470588235	20.588235294117645	5.88235294117647
[...]
