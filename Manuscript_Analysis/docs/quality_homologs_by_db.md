# Quality_homologs_by_db.py

Internal script
Run:

```
python3 uscripts/quality_homologs_by_db.py
```

### : Description :

Compute the proportion of homologs detected against each database passing 
the quality tests, comparing the proportions of total AMPs passing the same
tests in the whole AMPSphere.

### : Inputs :

1. **data/dramp_candidates.txt.gz**

List of AMP candidates matching to DRAMP v3 database with a significant E-value.
Based in MMSeqs2 results.

2. **data/gmgc_candidates.txt.gz**

List of AMP candidates matching to GMGC v1 database with a significant E-value.
These candidates were pre-filtered by using only AMP candidates passing the
terminal placement test - non-fragments.
Based in MMSeqs2 results.

3. **data/SmProt_candidates.txt.gz**

List of AMP candidates matching to SmProt2 database with significant E-value.
Based in MMSeqs2 results.

4. **data/starPepDB_candidates.txt.gz**

List of AMP candidates matching to starPep45k database with significant E-value.
Based in MMSeqs2 results.

5. **data/STsORFs_candidates.txt.gz**

List of AMP candidates matching to STsORFs database with significant E-value.
STsORFs are high confidence small ORFs from *Salmonella* genus.
Based in MMSeqs2 results.

6. **data/quality_assessment.tsv.gz**

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

1. **quality_by_db.tsv**

Table generated with the proportions of orthologs by database, 
passing in the different tests.

Columns represent the databases:
DRAMPv3, GMGCv1, SmProtv2, StarPepDB45k, and STsORFs.

Rows represent the quality tests:
Antifam, RNAcode, Metaproteomes, Metatranscriptomes, Terminal position

2. **quality_db_homologs.{svg, png}**

Vertical bar chart with `x axis` representing the quality tests, colors
representing the databases and `y axis` representing the proportion of 
homologs passed in the test.


