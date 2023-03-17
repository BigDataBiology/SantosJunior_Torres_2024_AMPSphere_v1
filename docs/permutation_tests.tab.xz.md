**permutation_tests.tab**

**Description:**	TSV-file containing the permutation results for protein families (full-length, all_fams) 
                        and less than 50 residues (short_fams). 50 thousand families were sampled randomly for
                        10 thousand times and the amount of families affiliating to conserved gene neighborhoods
                        were recorded. We also recorded those with genome contexts involving antibiotic synthesis
                        genes and antibiotic resistance (from CARD).

| **columns** | **description** |
| :---: | :---: |
| all_fams | families of proteins (all lengths) with conserved genome contexts | 
| all_fams_antibiotic | families of proteins (all lengths) with conserved genome contexts involving genes of antibiotic synthesis |
| all_fams_CARD | families of proteins (all lengths) with conserved genome contexts involving ARGs from CARD |
| short_fams | families of protein with a maximum of 50 residues with conserved genome contexts | 
| short_fams_antibiotic | families of protein with a maximum of 50 residues with conserved genome contexts involving genes of antibiotic synthesis |
| short_fams_CARD | families of protein with a maximum of 50 residues with conserved genome contexts involving ARGs from CARD |


**MD5 SUM:**	eaa18845442fbc10f39fe38ab78b1f99

**Size (MBytes):**	0.071596

**Content sample (first 5 items):**

all_fams	all_fams_antibiotic	all_fams_CARD	short_fams	short_fams_antibiotic	short_fams_CARD
18923	179	736	1964	30	124
19001	216	721	1930	31	119
19059	196	723	1935	24	118
18791	185	682	1835	28	99
19058	193	715	1961	35	105
18962	210	783	1934	30	116
18976	202	743	1971	28	134
18875	195	772	1970	17	118
19156	191	728	1965	28	121
[...]
