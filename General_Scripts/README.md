# AMPSphere_manuscript

Brings all code needed to regenerate AMPSphere, create the figures in the paper and perform the side-analysis.

|**Folder**|**Description**|
|:---:|:---:|
|AMPSphere_generation_v.2022-03|Folder used to regenerate AMPSphere from the base files|
|metadata_analysis|Folder used to aggregate metadata to the AMPSphere original files|
|figure_sup1|Folder used to generate the supplementary Figure 1 (Here it is performed the quality analysis and also the unique gene counting)|
|select_candidates_for_synthesis|Folder used to select candidates for synthesis, based in their solubility and in some synthesis empirical rules. It uses the quality of the candidates to curate a set of roughly 500 AMPs.|
|figure_1|Folder used to generate the Figure 1 with homology search results|
|figure_2|Folder used to generate the supplementary Figure 2, with cluster species affiliations|

Inside each folder the strucutre is an utils/ folder and a main.py script.
By executing this main.py script, automatically it will be generated two folders: analysis/ and data/
The data needed for the analysis will be downloaded automatically from the sources and the results will be available in the analysis/ folder.

Basically, to obtain the results, you should enter in the folder and execute:

```
    $ cd AMPSphere_generation_v.2022-03/
    $ python3 main.py
```

Detailed information about scripts, inputs and outputs are available in the README.md files inside each folder.

