**result_dramp.m8**

**Description:**	MMSeqs2 easy-search results using candidate AMPs from AMPSphere against the DRAMP v.3 database after quality filtering
                        leaving only mhits with identities above 75% covering at least 50% of the shorter sequence and with e-value < 1e-5.

| **column position** | **description** |
| :---: | :---: |
| 1 | AMP access code from the AMPSphere resource |
| 2 | AMP access code from the DRAMP v.3 database |
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

**MD5 SUM:**	2e7a983b7c530b71f2eeb82b2cf0616a

**Size (MBytes):**	1.6389694213867188

**Content sample (first 5 items):**

AMP10.273_700	DRAMP18659_GENERAL	2.907E-04	0	70.000	14	2	21	25	1	20	32	20	67	31	20M	AVRVAINGFGRIGRLAFRQMFGAEG	VKVGINGFGRIGRLVTRAAFHGKKVEVVAIND	AMP10.273_700 | SPHERE-III.004_469	DRAMP18659_GENERAL	VRVAINGFGRIGRLAFRQMF	VKVGINGFGRIGRLVTRAAF			6	0.800	0.625
AMP10.273_700	DRAMP18660_GENERAL	2.907E-04	0	70.000	14	2	21	25	1	20	32	20	67	31	20M	AVRVAINGFGRIGRLAFRQMFGAEG	VKVGINGFGRIGRLVTRAAFHGKKVEVVAIND	AMP10.273_700 | SPHERE-III.004_469	DRAMP18660_GENERAL	VRVAINGFGRIGRLAFRQMF	VKVGINGFGRIGRLVTRAAF			6	0.800	0.625
AMP10.292_996	DRAMP14850_PATENTED	8.169E-11	0	81.400	22	1	27	30	58	84	86	27	113	49	27M	ITPRGEKKAYVKLKPEYKASELAVKLGVIV	MDPYKVIIRPVVTEKAISLIEKENKLTFIVDRRATKTDIKKAIEEIFNVKVEKVNTLITPKGEKKAYVKLKPEYSASEIAARLGLF	AMP10.292_996 | SPHERE-III.016_390	DRAMP14850_PATENTED	ITPRGEKKAYVKLKPEYKASELAVKLG	ITPKGEKKAYVKLKPEYSASEIAARLG			5	0.900	0.314
AMP10.292_996	DRAMP14997_PATENTED	4.904E-08	0	67.800	19	1	28	30	58	85	86	28	95	42	28M	ITPRGEKKAYVKLKPEYKASELAVKLGVIV	MDAFDVIKAPVVTEKTVRMIEEENKLVFYVDRRATKQDIKRAMKELFDVEVEKVNTLITPKGEKKAYVKLKEGYDASKIAASLGIY	AMP10.292_996 | SPHERE-III.016_390	DRAMP14997_PATENTED	ITPRGEKKAYVKLKPEYKASELAVKLGV	ITPKGEKKAYVKLKEGYDASKIAASLGI			9	0.933	0.326
AMP10.617_764	DRAMP15338_PATENTED	3.641E-07	0	52.900	18	3	36	37	54	87	88	34	90	40	34M	ALAEPELMQAAQRNIIHKNNASRKVSRLAHQIAKLAK	MANTTSAKKATRKIARRSAVNKARRSRIRSFVRKVEEAIASGDQALAAAALKAAQPELMRPATKGVMHSNTASRKVSRLAQRVKSLSA	AMP10.617_764 | SPHERE-III.001_940	DRAMP15338_PATENTED	AEPELMQAAQRNIIHKNNASRKVSRLAHQIAKLA	AQPELMRPATKGVMHSNTASRKVSRLAQRVKSLS			16	0.919	0.386
[...]
