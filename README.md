# AMPSphere v2022-03

- This repository contains files and scripts to generate analysis and figures in AMPSphere
  paper.

## Info

The folder `General_Scripts` contains scripts needed to regenerate the AMPSphere resource from the
raw data. To have more info, follow the [link](General_Scripts/README.md).

The folder `Manuscript_Analysis` contains files needed to run the analysis and tests in the AMPSphere paper,
for more information, go to its [README.md](Manuscript_Analysis/README.md).

## Install

To install the conda environment needed to run these scritps, you can use conda:

```
  $ conda env create -f environment.yml -n ampsphere_manuscript
  $ conda activate ampsphere_manuscript
```

## Databases

These scripts rely on a series of pre-computed files and core resources that
are necessary to the proper functioning. Download them from [Zenodo](https://doi.org/10.5281/zenodo.7742544),
to this folder and then do following commands:

```
   $ tar -xf data_folder.tar
   $ tar -xf docs.tar
   $ rm -rf data_folder.tar docs.tar
```
