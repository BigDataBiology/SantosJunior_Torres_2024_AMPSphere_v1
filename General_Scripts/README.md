# AMPSphere_manuscript

Brings all code needed to regenerate AMPSphere resource:

[AMPSphere_generation_v.2022-03](AMPSphere_generation_v.2022-03/README)

    Folder used to regenerate AMPSphere from the base files
    
[metadata_analysis](metadata_analysis/README.md)

    Folder used to aggregate metadata to the AMPSphere original files

Inside each folder the strucutre is an `utils/` folder and a main.py script.

By executing this `main.py` script, automatically it will be generated two
folders: `analysis/` and `data/`.

The data needed for the analysis will be downloaded automatically from the
sources and the results will be available in the `analysis/` folder.

Detailed information about scripts, inputs and outputs are available in
the README.md files inside each folder.

