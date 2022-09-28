**ensemble_predictions.tsv**

**Description:**	During the process of selecting high-quality candidates for the *in vitro* testing,
                        the best candidates were submitted to different AMP prediction systems (AI4AMP, AmPEPpy, 
                        ampir, AMPlify, AMPScannerv2, APIN) to obtain those sequences co-predicted strengthening
                        our confidence. This TSV-file report the results for these systems and the final prediction
                        obtained for the candidates -- to that we adopted the thumb rule that >=50% of co-prediction
                        as an AMP would determine a true AMP.

| **columns** | **description** |
| :---: | :---: |
| amp | AMP candidate from AMPSphere -- access code |
| AI4AMP | result of prediction, either AMP or nonAMP |  # https://github.com/LinTzuTang/AI4AMP_predictor
| AmPEPpy | result of prediction, either AMP or nonAMP |  # https://github.com/tlawrence3/amPEPpy
| ampir | result of prediction, either AMP or nonAMP |  # https://ampir.marine-omics.net/
| AMPlify | result of prediction, either AMP or nonAMP |  # https://github.com/bcgsc/AMPlify
| AMPScannerv2 | result of prediction, either AMP or nonAMP |  # https://www.dveltri.com/ascan/
| APIN | result of prediction, either AMP or nonAMP |  # https://github.com/zhanglabNKU/APIN
| pcts | proportion of positive results |
| prediction | result of co-prediction, either AMP or nonAMP |

**MD5 SUM:**	9b6a44f41175dfca283841a4515c5dcf

**Size (MBytes):**	0.02428436279296875

**Content sample (first 5 items):**

amp	AI4AMP	AmPEPpy	ampir	AMPlify	AMPScannerv2	APIN	pcts	prediction
AMP10.000_033	nonAMP	nonAMP	AMP	nonAMP	AMP	AMP	50.0	AMP
AMP10.000_133	nonAMP	nonAMP	AMP	nonAMP	AMP	AMP	50.0	AMP
AMP10.000_176	nonAMP	nonAMP	AMP	nonAMP	AMP	AMP	50.0	AMP
AMP10.000_209	nonAMP	AMP	nonAMP	nonAMP	nonAMP	nonAMP	16.666666666666668	nonAMP
[...]
