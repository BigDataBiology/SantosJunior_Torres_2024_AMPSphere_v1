**db_motif.tsv**

**Description:**	Database used for the annotation of motifs in the AMP candidates from AMPSphere. 
                        The motifs were encoded as regex rules from the original references used:
                        - [Huan et al. (2020)](https://doi.org/10.3389/fmicb.2020.582779)
                        - [Ruhanen et al. (2014)](https://doi.org/10.3389/fmicb.2014.00004)
                        The database was structured as a TSV-file with 3 columns:

| **columns** | **description** |
| :---: | :---: |
| motif | ReGex rule to identify the motif |
| name | Name for a motif (it can be repeated) |
| description | Function of the motif |

**MD5 SUM:**	ca22c546f55cba512c021a0c327dd249

**Size (MBytes):**	0.00182342529296875

**Content sample (first 5 items):**

motif	name	description
LSGGQ	ABC_motif	ATP binding cassette transporter motif
^.K.{0,2}K	adherence	Signal peptide motif
^.K..S	adherence	Signal peptide motif
^.{1,2}K.{1,2}L	adherence	Signal peptide motif
[...]
