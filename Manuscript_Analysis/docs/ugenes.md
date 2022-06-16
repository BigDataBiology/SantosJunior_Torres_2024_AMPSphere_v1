# Ugenes.py

Internal script
Run:

```
python3 uscripts/ugenes.py
```

### : Description :

Plot the number of AMP candidates associated to 1, 2 or more genes. 

### : Inputs :

1. **data/AMPSphere_v.2022-03.fna.xz:**

Fasta file containing the nucleotide sequences for the AMP candidates in 
AMPSphere. The arrangement of this fasta file has headers for the sequences
following the standard:

```
>GMSC10.SMORF.000_000_000_949 | AMP10.000_004
    |_ GMSC accession                |_ AMP accession
```

### : Outputs :

1. **figure_unique_genes.svg**

Vertical bar chart with the number of gene variants per AMP in the `x axis` and
the number of AMPs in the `y axis`

