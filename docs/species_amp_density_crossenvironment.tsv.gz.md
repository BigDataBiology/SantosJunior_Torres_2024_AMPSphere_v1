**species_amp_density_per_sample.tsv.gz**

**Description:**	TSV-file containing the average AMP density and its devation for each species
                        happening in at least 10 samples per environment, and in at least 2 environments.
                        The difference of the medians were tested with Mann-Whitney's U and the significant
                        results after adjustment with Bonferroni were kept in the table.

| **columns** | **description** |
| :---: | :---: |
| species | species name |
| env1 | environment 1 tested |
| env1_samples | number of samples in which species appeared from environment 1 | 
| amp_density_avg_e1 | average amp density for species in samples from environment 1 |
| amp_density_std_e1 | standard deviation of amp density for species in samples from environment 1 |
| env2 | environment 2 tested |
| env2_samples | number of samples in which species appeared from environment 2 | 
| amp_density_avg_e2 | average amp density for species in samples from environment 2 |
| amp_density_std_e2 | standard deviation of amp density for species in samples from environment 2 |
| p_value | P-value from Mann-Whitney's U test between medians of two environments |
| p_value_adj | Adjusted P-value after Bonferroni test |

**MD5 SUM:**	956aaa8e298cf2c93b9c3116e3c182e9

**Size (MBytes):**	0.463566

**Content sample (first 5 items):**

species	env1	env1_samples	amp_density_avg_e1	amp_density_std_e1	env2	env2_samples	amp_density_avg_e2	amp_density_std_e2	p_valuep_adj
0-14-0-80-60-11 sp002779455	groundwater	27	1.0655357478748304	0.5593373058517308	marine	25	0.9871864607549803	0.5698100494745072	0.6209544971387669	1.0
0-14-0-80-60-11 sp002779455	groundwater	27	1.0655357478748304	0.5593373058517308	soil	155	1.2075633352987654	0.8342072291737697	0.7908466313996237	1.0
0-14-0-80-60-11 sp002779455	groundwater	27	1.0655357478748304	0.5593373058517308	water associated	45	0.8100397624649708	0.4050337955171445	0.06273530366791105	1.0
0-14-0-80-60-11 sp002779455	marine	25	0.9871864607549803	0.5698100494745072	soil	155	1.2075633352987654	0.8342072291737697	0.4466061500221108	1.0
0-14-0-80-60-11 sp002779455	marine	25	0.9871864607549803	0.5698100494745072	water associated	45	0.8100397624649708	0.4050337955171445	0.22967798475917311	1.0
0-14-0-80-60-11 sp002779455	soil	155	1.2075633352987654	0.8342072291737697	water associated	45	0.8100397624649708	0.4050337955171445	0.020818857476534995	1.0
0-14-3-00-41-53 sp002780895	groundwater	55	0.6977572406351122	0.2787958276378088	water associated	11	0.911448711134896	0.3808383435590022	0.05840564160941046	1.0
12-FULL-67-14b sp001768615	activated sludge	27	2.8808223094974146	1.423401597210508	plant associated	82	0.4194865711918924	0.2558606415723752	2.685626639382069e-14	2.5330830459443756e-10
12-FULL-67-14b sp001768615	activated sludge	27	2.8808223094974146	1.423401597210508	soil	213	1.007816662066981	1.1671434765936781	3.5342366465959486e-10	3.139457485618995e-06
[...]
