# Searching for Mycoplasma pneumoniae types

Use the jupyter notebook to download the genomes and them transfer those in the genomes_list to the folder 'genomes':

```
mkdir data/genomes
mv *fna.gz data/genomes/
```

After that, open the folder 'map' and run the pipeline as follow:

```
cd map/
bash pipeline_map.sh
```

Use the newick file generated: "p1_genes.nwk" at the [iTOL website](itol.embl.de/)
to visualize the groups and save their node ids.

Finally, after saving them into separate groups by strain type,
compare with the file "fam_red.tsv", which contains the
clusters per taxa.

