# Ftests.py

Internal script
Run:

```
python3 uscripts/ftests.py
```

### : Description :

It takes merged tables by samples and species of MAGs pre-processed using Macrel with AMP
candidates info. Then it screens out the MAGs generating cuts by AMPs per MAG and the outcome
and calculate the Exact Fisher's test for these partitions. It reports graphs and the 
tables contaning the odds ratio and the p-values. 

### : Inputs :

1. **data/merged_to_donor_samples.tsv.gz**

Merged information of preprocessed MAGs from Schmidt et al. (2022).
MAGs were pre-processed using Macrel and their information was merged to
the FMT.id, timepoint, sample (separately for donor, recipient pre- and post-FMT)
and species. Its columns are:

    species - cluster representing genome species
    fmt.id - id for the FMT
    fmt.type - 'FMT'
    timepoint.fmt_x - timepoint in the FMT procedure for the recipient sample
    sample.donor - biosample for the donor in this FMT procedure
    sample.pre - biosample for the recipient pre-FMT
    sample.post - biosample for the recipient post-FMT
    outcome - species outcome after FMT, e.g. retained, coexisted, lost
    sp_retained - float
    reject_consp - float
    coexist_consp.rec - float
    coexist_consp.donor - float
    engraft_consp - float
    reject_novel - float
    engraft_novel - float 
    sp_lost - float
    influx_consp - float
    influx_novel - float
    influx_total - float
    genome - genome name
    completeness - CheckM completeness
    contamination - CheckM contamination
    subject_type - origin (donor, recipient pre- or post-FMT)
    timepoint.fmt_y - timepoint in the FMT procedure for the genome sample 
    clinical_response - responder / non_responder
    ORFs - total number of ORFs (containing small ORFs) in the genome
    smORFs - total number of small ORFs predicted in the genome
    #_amps - number of non-redundant AMPs found in the genome
    sequences - non-redundant AMP sequences identified in the genome

2. **data/merged_to_recipient_samples_pre.tsv.gz**

Merged information of preprocessed MAGs from Schmidt et al. (2022).
MAGs were pre-processed using Macrel and their information was merged to
the FMT.id, timepoint, sample (separately for donor, recipient pre- and post-FMT)
and species. Its columns are:

    species - cluster representing genome species
    fmt.id - id for the FMT
    fmt.type - 'FMT'
    timepoint.fmt_x - timepoint in the FMT procedure for the recipient sample
    sample.donor - biosample for the donor in this FMT procedure
    sample.pre - biosample for the recipient pre-FMT
    sample.post - biosample for the recipient post-FMT
    outcome - species outcome after FMT, e.g. retained, coexisted, lost
    sp_retained - float
    reject_consp - float
    coexist_consp.rec - float
    coexist_consp.donor - float
    engraft_consp - float
    reject_novel - float
    engraft_novel - float 
    sp_lost - float
    influx_consp - float
    influx_novel - float
    influx_total - float
    genome - genome name
    completeness - CheckM completeness
    contamination - CheckM contamination
    subject_type - origin (donor, recipient pre- or post-FMT)
    timepoint.fmt_y - timepoint in the FMT procedure for the genome sample 
    clinical_response - responder / non_responder
    ORFs - total number of ORFs (containing small ORFs) in the genome
    smORFs - total number of small ORFs predicted in the genome
    #_amps - number of non-redundant AMPs found in the genome
    sequences - non-redundant AMP sequences identified in the genome

3. **data/merged_to_recipient_samples_post.tsv.gz**

Merged information of preprocessed MAGs from Schmidt et al. (2022).
MAGs were pre-processed using Macrel and their information was merged to
the FMT.id, timepoint, sample (separately for donor, recipient pre- and post-FMT)
and species. Its columns are:

    species - cluster representing genome species
    fmt.id - id for the FMT
    fmt.type - 'FMT'
    timepoint.fmt_x - timepoint in the FMT procedure for the recipient sample
    sample.donor - biosample for the donor in this FMT procedure
    sample.pre - biosample for the recipient pre-FMT
    sample.post - biosample for the recipient post-FMT
    outcome - species outcome after FMT, e.g. retained, coexisted, lost
    sp_retained - float
    reject_consp - float
    coexist_consp.rec - float
    coexist_consp.donor - float
    engraft_consp - float
    reject_novel - float
    engraft_novel - float 
    sp_lost - float
    influx_consp - float
    influx_novel - float
    influx_total - float
    genome - genome name
    completeness - CheckM completeness
    contamination - CheckM contamination
    subject_type - origin (donor, recipient pre- or post-FMT)
    timepoint.fmt_y - timepoint in the FMT procedure for the genome sample 
    clinical_response - responder / non_responder
    ORFs - total number of ORFs (containing small ORFs) in the genome
    smORFs - total number of small ORFs predicted in the genome
    #_amps - number of non-redundant AMPs found in the genome
    sequences - non-redundant AMP sequences identified in the genome
    
### : Outputs :

1. **fisher_test_{cutoff}.tsv**

Using different cutoffs for the number of AMP candidates per MAG, the
table shows the group, outcome and the Exact Fisher's test for the
partitions performed. It returns the outcomes significantly influenced by the
number of AMPs carried in a MAG. Its columns are:

    group - donor, recipient pre- or post-FMT
    outcome - species outcome after FMT procedure
    odds_ratio - float, proportion of the two groups ('<x', '>=x')
    pvalue - float, probability 

2. **fisher_test_summary.tsv**

Table containing all results from the tests using different cutoffs.
Its columns are:

    group - donor, recipient pre- or post-FMT
    outcome - species outcome after FMT procedure
    odds_ratio - float, proportion of the two groups ('<x', '>=x')
    pvalue - float, probability
    amp_cutoff - number of AMPs per MAG to perform the data split

3. **donor_{cutoff}_test.svg**

Related to the merged samples of donor. Bar chart with the proportion of an
outcome in the `y axis` and the outcome in the `x axis`. Bars related to the
standard error are shown, and were calculated as: np.sqrt(p(1-p)/n)
Colors represent groups accordingly the cutoff (AMPs per MAG) established:
'<x', '>=x'. Only outcomes significantly affected by the AMPs were shown.

4. **pre-FMT_{cutoff}_test.svg**

Related to the merged samples of recipient pre-FMT. Bar chart with the proportion
of an outcome in the `y axis` and the outcome in the `x axis`. Bars related to the
standard error are shown, and were calculated as: np.sqrt(p(1-p)/n)
Colors represent groups accordingly the cutoff (AMPs per MAG) established:
'<x', '>=x'. Only outcomes significantly affected by the AMPs were shown.

5. **post-FMT_{cutoff}_test.svg**

Related to the merged samples of recipient post-FMT. Bar chart with the proportion
of an outcome in the `y axis` and the outcome in the `x axis`. Bars related to the
standard error are shown, and were calculated as: np.sqrt(p(1-p)/n)
Colors represent groups accordingly the cutoff (AMPs per MAG) established:
'<x', '>=x'. Only outcomes significantly affected by the AMPs were shown.

