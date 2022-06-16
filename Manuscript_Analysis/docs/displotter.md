# Displotter.py

Internal script
Run:

```
python3 uscripts/displotter.py
```

### : Description :

It generates a graph with the data for the validation of the clustering
procedures.

### : Inputs :

1. **data/output_clustering_significance_levelI.tsv.gz:**

Table from AMPSphere v2022-03 version. It consists in a table generated during
the validation of the clustering procedure. The validation consisted in the
alignment of 1000 randomly selected AMP candidates against the respective
representatives of their clusters in a given level.

The alignment was performed using Smith-Waterman algorithm in triplicate and
the identity/significance were calculated as in blast. Its columns are:
    
    query - randomly AMP candidate (excluded cluster representatives)
    target - cluster representative
    identity - identity calculated counting gaps
    gap_identity - identity calculated excluding gaps
    coverage - query coverage in percent of residues
    score - alignment score
    aln_len - alignment length
    evalue - expected-value from the exponential distribution by Altschull et al.
    sig. - '*' if e-value < 1e-5 else 'n.s.'
    family - AMP cluster accession at 1st level (100% identity)
    replicate - replicate batch (1, 2, or 3)


2. **data/output_clustering_significance_levelII.tsv.gz**

Table from AMPSphere v2022-03 version. It consists in a table generated during
the validation of the clustering procedure. The validation consisted in the
alignment of 1000 randomly selected AMP candidates against the respective
representatives of their clusters.

The alignment was performed using Smith-Waterman algorithm in triplicate and
the identity/significance were calculated as in blast. Its columns are:
    
    query - randomly AMP candidate (excluded cluster representatives)
    target - cluster representative
    identity - identity calculated counting gaps
    gap_identity - identity calculated excluding gaps
    coverage - query coverage in percent of residues
    score - alignment score
    aln_len - alignment length
    evalue - expected-value from the exponential distribution by Altschull et al.
    sig. - '*' if e-value < 1e-5 else 'n.s.'
    family - AMP cluster accession at 2nd level (100-85% identity)
    replicate - replicate batch (1, 2, or 3)

3. **data/output_clustering_significance_levelIII.tsv.gz**

Table from AMPSphere v2022-03 version. It consists in a table generated during
the validation of the clustering procedure. The validation consisted in the
alignment of 1000 randomly selected AMP candidates against the respective
representatives of their clusters.

The alignment was performed using Smith-Waterman algorithm in triplicate and
the identity/significance were calculated as in blast. Its columns are:
    
    query - randomly AMP candidate (excluded cluster representatives)
    target - cluster representative
    identity - identity calculated counting gaps
    gap_identity - identity calculated excluding gaps
    coverage - query coverage in percent of residues
    score - alignment score
    aln_len - alignment length
    evalue - expected-value from the exponential distribution by Altschull et al.
    sig. - '*' if e-value < 1e-5 else 'n.s.'
    family - AMP cluster accession at 3rd level (100-85-75% identity)
    replicate - replicate batch (1, 2, or 3)

### : Outputs :

1. **clustering_significance.svg**

3 scatterplots disposed horizontally. They represent respectively the levels I, II, and III of
clustering used in the AMPSphere. In their `x axis` is the identity and `y axis` is the Log(E-value)
for the alignment of a randomly selected protein against the representative of the cluster it belongs
to. The color of the dots represent a replicate of the experiment containing 1000 AMP candidates.

