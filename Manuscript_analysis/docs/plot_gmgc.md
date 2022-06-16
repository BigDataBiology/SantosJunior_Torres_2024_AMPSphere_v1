# Plot_gmgc.py

Internal script
Run:

```
python3 uscripts/plot_gmgc.py
```

### : Description :

It plots the results for the homologs search using AMPSphere against GMGC v1.

### : Inputs :

1. **data/result_gmgc.m8.xz:**

Table of homologs search using AMPSphere v2022-03 version against GMGC v1.
Only candidates passing the terminal placement test, in other words, 
not fragments. It was generated by MMSeqs2 and consists of:
    
    query - amp accession in AMPSphere v2022-03 
    target - GMGC target sequence accession
    evalue - expected-value
    gapopen - number of gap open
    pident - protein identity
    nident - N identity
    qstart - alignment start in the query
    qend - alignment end in the query
    qlen - length of query
    tstart - alignment start in the target
    tend - alignment end in the target
    tlen - length of target
    alnlen - alignment length
    raw - x
    bits - bit score
    cigar - CIGAR hash
    qseq - query sequence
    tseq - target sequence
    qheader - description of query AMP
    theader - description of target GMGC sequence
    qaln - aligned query sequence
    taln - aligned target sequence
    qframe - x
    tframe - x
    mismatch - number of mismatches 
    qcov - query coverage
    tcov - target coverage

### : Outputs :

1. **out_test.svg**

Figure with 4 panes:

    a) Scatterplot, x axis = identity, y axis = Log(E-value)
    b) Scatterplot, x axis = length of AMP, y axis = GMGC protein length
    c) KDE plot, x axis = start site of alignment in GMGC protein (as %), hue = Identity ('*')
    d) KDE plot, x axis = Query coverage (%), hue = Identity ('*')
    e) KDE plot, x axis = Target coverage (%), hue = Identity ('*')
    f) Bar chart, x axis = GMGC ortholog group, y axis = Log(counts)
            
'*' Identity classes: 0-25%, 25%-50%, 50%-75%, 75%-100%
