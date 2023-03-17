# Transmissibility and AMP density

The transmissibility uses the species with a minimum prevalence of
10 samples per environment and cross data with the transmissibility
measured in species detected from gut and mouth in humans between
mother and offspring at different ages that you can find in this
[link](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9892008/).

To reproduce the data used in this analysis:

```
    $ python3 scripts/test_species_content.py
```

Detailed inputs list:

data = pd.read_table('data/')
metadata = pd.read_table("data/")


| **Input file** | **Description** |
| :---: | :---: |
| complete_amps_associated_taxonomy.tsv.gz | table generated during metadata analysis with columns amp, gene, sample, source (species), is_metagenomic (bool: True/False) |
| metadata.tsv.xz | contains the data for the samples, including their sampling material and location |

Detailed outputs list:

test.to_csv('analysis/',

| **Input file** | **Description** |
| :---: | :---: |
| species_amp_richness_crossenvironment.tsv.gz | detailed information about the FMTs and their outcomes with the MAGs species per event |


