## AMPSphere v.2022-03 [core analysis]

This folder brings the different steps to generate the files
used in the analysis made and shown in the AMPSphere resource
paper. Among the folders there are those needed to generate
the data in the AMPSphere as well as those needed to rearrange
data needed for those steps, such as working the metadata
and merging information.

### Data

All data needed throughout the analysis are linked as 
soft links in the respective folders. The structure for each
analysis is basically in three folders:

1. data - contains the links to the data used in the analysis
2. utils - contains the scripts for the different analysis
3. analysis - where the results will be saved

It is also available a **README.md** file with general and 
specific information of the analysis goals, steps, input and
output files as well as brief info for the scripts and functions.

In a general way, the script saved in the analysis folder as 
*main.py* is the only command needed to execute the entire action
by typing:

```
$ python main.py
```
The files used in all analysis are saved in the `data_folder`,
which the contents are listed in the link below:

[data_folder](data_folder/README.md)

### Analysis

To verify each analysis in this folder, you can quickly 
follow the links below to move to their functionalities 
and specificities.

 - [01_resource_generation](01_resource_generation/README.md)
 - [02_motifs_annotation](02_motifs_annotation/README.md)
 - [03_clustering_significance](03_clustering_significance/README.md)
 - [04_metadata_analysis](04_metadata_analysis/README.md)
 - [05_clonality_and_quality](05_clonality_and_quality/README.md)
 - [06_rnacode_amp_families](06_rnacode_amp_families/README.md)
 - [07_homologs_and_overlap](07_homologs_and_overlap/README.md)
 - [08_BGI_peptides](08_BGI_peptides/README.md)
 - [09_taxonomy_core_analysis](09_taxonomy_core_analysis/README.md)
 - [10_select_candidates_for_synthesis](10_select_candidates_for_synthesis/README.md)
 - [11_AMPs_from_Mpneumoniae](11_AMPs_from_Mpneumoniae/README.md)
 - [12_host_non_host_amps](12_host_non_host_amps/README.md)
 - [13_fmt_story](13_fmt_story/README.md)

A complete relationship of input and output files
per analysis is available [here](docs/Input_vs_Ouput_per_analysis.md)

