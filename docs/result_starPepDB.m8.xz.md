**result_starPepDB.m8**

**Description:**	MMSeqs2 easy-search results using candidate AMPs from AMPSphere against
                        the [StarPepDB 45k database](http://mobiosd-hub.com/starpep/) after quality
                        filtering leaving only mhits with identities above 75% covering at least
                        50% of the shorter sequence and with e-value < 1e-5.

| **column position** | **description** |
| :---: | :---: |
| 1 | AMP access code from the AMPSphere resource |
| 2 | AMP access code from the StarPepDB 45k database |
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


**MD5 SUM:**	14784018d7404213f3d6729dfee1291d

**Size (MBytes):**	1.6447372436523438

**Content sample (first 5 items):**

AMP10.161_382	starPep_30262	2.628E-07	0	65.500	19	1	29	34	52	80	98	29	93	41	29M	GLGKDYTLFALIDGYVKFESFKGKKRVSIYAEKK	MASKASGGSTRNGRDSISKRLGVKRYDGQFVKAGNIIVRQRGTRIYPGKNVGMGSDYTLFALKDGYVYFETRRKKKFVSVLSPEEWEKVMAQKNGKVH	AMP10.161_382 | SPHERE-III.012_440	starPep_30262	GLGKDYTLFALIDGYVKFESFKGKKRVSI	GMGSDYTLFALKDGYVYFETRRKKKFVSV			10	0.853	0.296
AMP10.161_382	starPep_29903	6.880E-07	1	54.500	18	1	32	34	52	84	88	33	90	40	19M1D13M	GLGKDYTLFALIDGYVKFESFKGKKRVSIYAEKK	MAHKKGQGSTQNNRDSAGRRLGVKKFGSEFVRAGNIIVRQRGTKMHPGNNVGMGKDHTLYALIDGVVKFEHKDRNRKKVSVVSQNFGE	AMP10.161_382 | SPHERE-III.012_440	starPep_29903	GLGKDYTLFALIDGYVKFE-SFKGKKRVSIYAE	GMGKDHTLYALIDGVVKFEHKDRNRKKVSVVSQ			14	0.941	0.375
AMP10.161_382	starPep_29905	1.802E-06	1	58.800	20	1	33	34	52	85	86	34	85	38	19M1D14M	GLGKDYTLFALIDGYVKFESFKGKKRVSIYAEKK	MAHKKGSGSTRNGRDSNSKRLGVKKYGGEQVTAGNILIRQRGTKVKPGQNVGKGKDDTLFALIDGFVLFEKSNQKQKTISVYSSKN	AMP10.161_382 | SPHERE-III.012_440	starPep_29905	GLGKDYTLFALIDGYVKFE-SFKGKKRVSIYAEK	GKGKDDTLFALIDGFVLFEKSNQKQKTISVYSSK			13	0.971	0.395
AMP10.161_382	starPep_29906	1.802E-06	1	61.200	19	1	30	34	52	82	87	31	85	38	19M1D11M	GLGKDYTLFALIDGYVKFESFKGKKRVSIYAEKK	MAHKKGTGSTRNGRDSNAQRLGVKRYGGQTVTAGSIIVRQRGTQVHPGNNVGRGKDDTLFALIDGVVKFEHKTRSRRKVSVYPATAE	AMP10.161_382 | SPHERE-III.012_440	starPep_29906	GLGKDYTLFALIDGYVKFE-SFKGKKRVSIY	GRGKDDTLFALIDGVVKFEHKTRSRRKVSVY			11	0.882	0.356
AMP10.161_382	starPep_32626	8.983E-06	1	61.200	19	1	30	34	60	90	94	31	80	36	21M1D9M	GLGKDYTLFALIDGYVKFESFKGKKRVSIYAEKK	MLRLDLQFFASKKGVGSTKNGRDSEAKRLGAKRADGQFVTGGSILYRQRGTKIYPGENVGRGGDDTLFAKIDGTVKFERFGRDRKKVSVYPVAQ	AMP10.161_382 | SPHERE-III.012_440	starPep_32626	GLGKDYTLFALIDGYVKFESF-KGKKRVSIY	GRGGDDTLFAKIDGTVKFERFGRDRKKVSVY			11	0.882	0.330
[...]
