## Non_host- and Host-associated AMP contents

This script set will create diverse tables related to the rarity of candidate AMPs 
in the different environments as well as the number of samples in which each
AMP occurs.

The inputs used are all described bellow:

| **Input file** | **Description** |
| :---: | :---: |
| uniques.tar.gz | this folder needs to be uncompressed before use, it consists of NGLess outputs of the mapping of genes from AMPSphere against the samples, these genes were only clustered at 100% identity and length prior analysis |
| envo.pkl | pickled dictionary containing two other dictionaries: 'general_envo_name' and 'high'. Each dictionary contains the sample as a key and the environment as a value | 
| headers.tsv.xz | table generated in the AMPSphere resource, it consists of the gene access code and the corresponding AMP |

The outputs are generated in the `analysis/` folder, and consist of:

| **Output file** | **Description** |
| :---: | :---: |
| sparse_matrix.npz | sparse matrix containing presence/absence as bool for AMPs in the rows and samples in the columns |
| saved_sample_cols.txt | key to decode columns in the sparse matrix |
| unique_map_cAMPs_high_habs.tsv.xz | output table containing the number of samples per high-level habitat in which each AMP was spotted |
| unique_map_cAMPs_general_habs.tsv.xz | output table containing the number of samples per low-level habitat in which each AMP was spotted |
| unique_map_cAMPs_nsamples.tsv.xz | output table containing the total number of samples in which each AMP was spotted |

The scripts used in this task are in the `utils/` folder and are explained bellow:

| **Script** | **Function** | **Input** | **Outputs** |
| :---: | :---: | :---: | :---: |
| make_all_one_28112022.py | generate tables of number of samples per habitat in which we find each AMPs from mapping data | All inputs | Outputs to anaysis/ all the outputs |

All scripts and functions are united in the main.py script and to generate the results
of this analysis, the user should simply:

```
    $ ./main.sh
```

