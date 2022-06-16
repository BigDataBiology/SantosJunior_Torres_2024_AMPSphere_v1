# Jelena_s.py

Internal script
Run:

```
python3 uscripts/jelena_s.py
```

### : Description :

It generates kde plots for different features of AMPSphere, Macrel traning set,
and DRAMP v3.0. 

### : Inputs :

1. **data/dramp_v3.features.tsv.gz:**

The AMPs from [DRAMP v3](DRAMP.cpu-bioinfor.org/) were filtered to eliminate the peptides
containing non-cannonical residues. 

This table contains features calculated with [Macrel](macrel.readthedocs.io/) pipeline.
Its columns are:
    
    Access - AMP accession in DRAMP v3
    sequence - amino acids sequence
    group - decoy column, label 'Unk'
    tinyAA - proportion of (A+C+G+S+T) in the peptide
    smallAA - proportion of (A+B+C+D+G+N+P+S+T+V) in the peptide
    aliphaticAA - proportion of (A+I+L+V) in the peptide
    aromaticAA - proportion of (F+H+W+Y) in the peptide
    nonpolarAA - proportion of (A+C+F+G+I+L+M+P+V+W+Y) in the peptide
    polarAA - proportion of (D+E+H+K+N+Q+R+S+T+Z) in the peptide
    chargedAA - proportion of (B+D+E+H+K+R+Z) in the peptide
    basicAA - proportion of (H+K+R) in the peptide
    acidicAA - proportion of (B+D+E+Z) in the peptide
    charge - calculated charge at pH 7 with pKscale "EMBOSS"
    pI - isoelectric point using pKscale "EMBOSS"
    aindex - relative volume occupied by aliphatic side chains (AVIL)
    instaindex - stability of a protein based on its amino acid composition
    boman - overall estimate of the potential of a peptide to bind to membranes or other
            proteins as receptor
    hydrophobicity - average hydrophobicity using scale "KyteDoolittle"
    hmoment - quantitative measure of the amphiphilicity perpendicular to the axis of any
              periodic peptide structure, such as the alpha-helix or beta-sheet
    SA.Group1.residue0 - CTDD at first residue, using solvent accessibility, class I
    SA.Group2.residue0 - CTDD at first residue, using solvent accessibility, class II
    SA.Group3.residue0 - CTDD at first residue, using solvent accessibility, class III
    HB.Group1.residue0 - CTDD at first residue, using FET from water to lipid phase, CI
    HB.Group2.residue0 - CTDD at first residue, using FET from water to lipid phase, CII
    HB.Group3.residue0 - CTDD at first residue, using FET from water to lipid phase, CIII

Solvent Accessibility (SA) was adopted as in previous studies (Dubchak et al. 1995, 1999),
however, the new feature HB was adapted from Von Heijne and Blomberg, 1979. Classes adopted
to the sequence encoding of the distribution at the first residue of each class:

 | **Properties** | **Class I** | **Class II** | **Class III** |
 | :---: | :---: | :---: | :---: |
 | SA | A, L, F, C, G, I, V, W | R, K, Q, E, N, D | M, S, P, T, H, Y |
 | FT | I,L,V,W,A,M,G,T | F,Y,S,Q,C,N | P,H,K,E,D,R |

2. **data/macrel_trainpos.features.tsv.gz:**

The AMPs from [Macrel traning set](https://github.com/BigDataBiology/macrel/tree/master/train)
were used to calculate features with [Macrel](macrel.readthedocs.io/) pipeline.
Its columns are:
    
    Access - AMP accession in Macrel training set
    sequence - amino acids sequence
    group - decoy column, label 'Unk'
    tinyAA - proportion of (A+C+G+S+T) in the peptide
    smallAA - proportion of (A+B+C+D+G+N+P+S+T+V) in the peptide
    aliphaticAA - proportion of (A+I+L+V) in the peptide
    aromaticAA - proportion of (F+H+W+Y) in the peptide
    nonpolarAA - proportion of (A+C+F+G+I+L+M+P+V+W+Y) in the peptide
    polarAA - proportion of (D+E+H+K+N+Q+R+S+T+Z) in the peptide
    chargedAA - proportion of (B+D+E+H+K+R+Z) in the peptide
    basicAA - proportion of (H+K+R) in the peptide
    acidicAA - proportion of (B+D+E+Z) in the peptide
    charge - calculated charge at pH 7 with pKscale "EMBOSS"
    pI - isoelectric point using pKscale "EMBOSS"
    aindex - relative volume occupied by aliphatic side chains (AVIL)
    instaindex - stability of a protein based on its amino acid composition
    boman - overall estimate of the potential of a peptide to bind to membranes or other
            proteins as receptor
    hydrophobicity - average hydrophobicity using scale "KyteDoolittle"
    hmoment - quantitative measure of the amphiphilicity perpendicular to the axis of any
              periodic peptide structure, such as the alpha-helix or beta-sheet
    SA.Group1.residue0 - CTDD at first residue, using solvent accessibility, class I
    SA.Group2.residue0 - CTDD at first residue, using solvent accessibility, class II
    SA.Group3.residue0 - CTDD at first residue, using solvent accessibility, class III
    HB.Group1.residue0 - CTDD at first residue, using FET from water to lipid phase, CI
    HB.Group2.residue0 - CTDD at first residue, using FET from water to lipid phase, CII
    HB.Group3.residue0 - CTDD at first residue, using FET from water to lipid phase, CIII

Solvent Accessibility (SA) was adopted as in previous studies (Dubchak et al. 1995, 1999),
however, the new feature HB was adapted from Von Heijne and Blomberg, 1979. Classes adopted
to the sequence encoding of the distribution at the first residue of each class:

 | **Properties** | **Class I** | **Class II** | **Class III** |
 | :---: | :---: | :---: | :---: |
 | SA | A, L, F, C, G, I, V, W | R, K, Q, E, N, D | M, S, P, T, H, Y |
 | FT | I,L,V,W,A,M,G,T | F,Y,S,Q,C,N | P,H,K,E,D,R |

3. **data/ampsphere_v2022-03.features.tsv.gz:**

The AMPs from [AMPSphere v2022-03](ampsphere.big-data-biology.org/)
were used to calculate features with [Macrel](macrel.readthedocs.io/) pipeline.
Its columns are:
    
    Access - AMP accession in Macrel training set
    sequence - amino acids sequence
    group - decoy column, label 'Unk'
    tinyAA - proportion of (A+C+G+S+T) in the peptide
    smallAA - proportion of (A+B+C+D+G+N+P+S+T+V) in the peptide
    aliphaticAA - proportion of (A+I+L+V) in the peptide
    aromaticAA - proportion of (F+H+W+Y) in the peptide
    nonpolarAA - proportion of (A+C+F+G+I+L+M+P+V+W+Y) in the peptide
    polarAA - proportion of (D+E+H+K+N+Q+R+S+T+Z) in the peptide
    chargedAA - proportion of (B+D+E+H+K+R+Z) in the peptide
    basicAA - proportion of (H+K+R) in the peptide
    acidicAA - proportion of (B+D+E+Z) in the peptide
    charge - calculated charge at pH 7 with pKscale "EMBOSS"
    pI - isoelectric point using pKscale "EMBOSS"
    aindex - relative volume occupied by aliphatic side chains (AVIL)
    instaindex - stability of a protein based on its amino acid composition
    boman - overall estimate of the potential of a peptide to bind to membranes or other
            proteins as receptor
    hydrophobicity - average hydrophobicity using scale "KyteDoolittle"
    hmoment - quantitative measure of the amphiphilicity perpendicular to the axis of any
              periodic peptide structure, such as the alpha-helix or beta-sheet
    SA.Group1.residue0 - CTDD at first residue, using solvent accessibility, class I
    SA.Group2.residue0 - CTDD at first residue, using solvent accessibility, class II
    SA.Group3.residue0 - CTDD at first residue, using solvent accessibility, class III
    HB.Group1.residue0 - CTDD at first residue, using FET from water to lipid phase, CI
    HB.Group2.residue0 - CTDD at first residue, using FET from water to lipid phase, CII
    HB.Group3.residue0 - CTDD at first residue, using FET from water to lipid phase, CIII

Solvent Accessibility (SA) was adopted as in previous studies (Dubchak et al. 1995, 1999),
however, the new feature HB was adapted from Von Heijne and Blomberg, 1979. Classes adopted
to the sequence encoding of the distribution at the first residue of each class:

 | **Properties** | **Class I** | **Class II** | **Class III** |
 | :---: | :---: | :---: | :---: |
 | SA | A, L, F, C, G, I, V, W | R, K, Q, E, N, D | M, S, P, T, H, Y |
 | FT | I,L,V,W,A,M,G,T | F,Y,S,Q,C,N | P,H,K,E,D,R |

### : Outputs :

1. **all_features.svg**

KDE plots showing the distribution of different features from the
AMPSphere, Macrel traning set and DRAMP. The panels follow below: 

 |  **label**  |  **feature**  |
 |  :---:  |  :---:  |
 |  a)  | # length (with legend)
 |  b)  | # smallAA  |
 |  c)  | # basicAA  |
 |  d)  | # pI  |
 |  e)  | # charge  |
 |  f)  | # aindex  |
 |  g)  | # instaindex  |
 |  h)  | # boman  |
 |  i)  | # hmoment  |

