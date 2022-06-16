# Qualtest.py

Internal script
Run:

```
python3 uscripts/qualtest.py
```

### : Description :

Plot the proportion of AMPs passing, failing or not tested
by quality test applied to AMPSphere.

### : Inputs :

1. **data/quality_assessment.tsv.gz**

Table containing the quality assessment per AMP in AMPSphere.
Its columns are:

    AMP - AMP accession code in AMPSphere
    Antifam - Match to Antifam (Failed) or not (Passed)
    RNAcode - Passed RNAcode test
    metaproteomes - Match to peptides from metaproteomes
    metatranscriptomes - Match to transcripts from at least 2
                         metatranscriptomes samples
    Coordinates - Is preceeded by a stop codon (non-fragment, Passed)
                  or not (potential fragment, Failed)


### : Outputs :

1. **figure_S1a_amp_quality.svg**

Horizontal bar chart with `x axis` with the number of AMP candidates
and `y axis` with the quality tests applied to the dataset.

The bars show the proportion of candidates passing, failing or
not tested. 

