# Homolog_quality.py

Internal script
Run:

```
python3 uscripts/homolog_quality.py
```

### : Description :

It plots the results of the mapping of AMP candidates against large proteins from GMGC
returning the number of mapped AMPs per database and their quality.


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

6. **data/quality_candidates.txt.gz:**

List of candidates passing simultaneously in all quality tests (RNAcode, AntiFam, and
Terminal Placement), but without experimental evidence of translation or transcription.
This list is generated during the quality assessment and can be reobtained by using
the file `quality_assessment.tsv.gz`.

7. **data/high_quality_candidates.txt.gz:**

List of candidates passing simultaneously in all quality tests (RNAcode, AntiFam, and
Terminal Placement) with experimental evidence of translation or transcription.
This list is generated during the quality assessment and can be reobtained by using
the file `quality_assessment.tsv.gz`.


### : Outputs :

1. **bar_chart_homologs_n_qual.svg**

Vertical bar chart with the number of AMPs homologs (in the `y axis`) to the databases
used to search (`x axis`). Gray codes for AMP candidates not passing in all tests and
black means the candidates that passed all quality control. The dashed red line shows
the expected number of high-quality candidates for a sample of the same size of candidates
randomly taken from AMPSphere.

2. **homologs_1.svg**

Pie chart with the proportion of AMP candidates matching to known databases and those
that remained unknown.


