## Non_host- and Host-associated AMP contents

This scripts set will create diverse tables related to the species from which AMPs 
are comming from in different environments along with other informations about the
contents of AMPs in host- and non-host-associated samples. The main output is a box-
plot with the different environments and their normalized AMP contents.

The inputs used are all described bellow:

| **Input file** | **Description** |
| :---: | :---: |
| gmsc_amp_genes_envohr_source.tsv.gz | table generated in metadata analysis, consists of the: sample, source, specI, is_metagenomic, geographic_location, latitude, longitude, general_envo_name, environment_material associated to genes and amps |
| samples-min500k-assembly-prodigal-stats.tsv | table generated in the GMSC resource, it consists of information about the assemblies and gene prediction per sample | 
| SPHERE_v.2022-03.levels_assessment.tsv.gz | table generated in the AMPSphere resource, it consists of the AMPs and their clustering information |
| taxonomy_annotation.tsv | table generated during the analysis for the Figure 2, it consists of: amp, taxid, level, source, fixed |

The outputs are generated in the `analysis/` folder, and consist of:

| **Output file** | **Description** |
| :---: | :---: |
| all_amps_hostvsenv.tsv | output bringing the information about the AMPs per sample, their aseembly length and the normalized value |
| shannon_diversity_per_sample.tsv | output table bringing the shannon diversity for the amp and the species as microbial sources |
| amps_list_allsamples.txt | amp lists for each macroenvironment tested and can be used as input for the Interactivenn software |
| species_distribution_by_env.tsv | output table containing the species and the number of AMPs per different tested macroenvironments |
| normamps_offoutlier_hostvsenv.svg | boxplot bringing the normalized number of AMPs per assembly length per different tested macroenvironments |

The scripts used in this task are in the `utils/` folder and are explained bellow:

| **Script** | **Function** | **Input** | **Outputs** |
| :---: | :---: | :---: | :---: |
| download_files.py | Downloads to folder `data/`the inputs needed for the main.py execution | None | Outputs to data/ all the inputs |
| infolists.py | Contains the lists of microontology terms per macroenvironment tested | None | None |
| general_functions.py | Contains functions to diverse uses for different scripts | None | None |
| load_data.py | Loads and pre-treat the inputs | gmsc_amp_genes_envohr_source.tsv.gz, SPHERE_v.2022-03.levels_assessment.tsv.gz, taxonomy_annotation.tsv | 2 dataframes (data, spheres) |
| distributions.py | Calculates distribution and diversity of species and AMPs per sample | 2 dataframes (data, spheres) | shannon_diversity_per_sample.tsv, species_distribution_by_env.tsv |
| ampsets_generator.py | Generate sets of AMPs per sample and count them for each macroenvironment | dataframe(data) | exports the files all_amps_hostvsenv.tsv, amps_list_allsamples.txt, and generates the dataframe (df) |
| boxstats.py | Generates the boxplot for normalized measures of AMP/sample and test their difference with MannWhitney-U test | dataframe(df) | normamps_offoutlier_hostvsenv.svg |

All scripts and functions are united in the main.py script and to generate the results
of this analysis, the user should simply:

```
    $ python3 main.py
```

