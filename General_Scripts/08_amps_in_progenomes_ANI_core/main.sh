#!/bin/bash

# to generate reference AMPs for progenomes 
python3 utils/generating_AMP_progenomes.py

# to download genomes (nucleotides)
python3 utils/ani_brutte.py  

# to cluster clones and strains
# to analyze the statistics behind AMPs shared by strains
python3 utils/getting_clusters.py

# to download proteomes 
python3 utils/protein_brutte.py

# to cluster and test prevalence of full-length protein families
python3 utils/workfams.py

# to analyze the core protein families
python3 utils/analysis_core_fams_res.py
