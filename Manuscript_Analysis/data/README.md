# AMPSphere v2022-03 - Analysis

## Data for manuscript resources

These files are the inputs for the analysis in the paper for
AMPSphere v2022-03. Their detailed description follow below.

### Files

1. **data/pgenomes_samples.tsv.gz:**

Resource from ProGenomes2 database, it consists in the
samples and their respective specI clusters, e.g.:

```
specI_v3_Cluster1	1049759.SAMN02436435
specI_v3_Cluster2	1001589.SAMN02436335
```

2. **data/gmsc_amp_genes_envohr_source.tsv.gz:**

Table from AMPSphere v2022-03 version. It consists in a table linking the
GMSC10 genes and the AMPs in the final resource. This table also contains
general metadata information. Its columns are:
    
    gmsc - gene accession in GMSC10 resource
    amp - AMP accession in AMPSphere v2022-03 resource
    sample - biosample from which we could identify the gene
    source - contig taxonomy from GTDB or ProGenomes2
    specI - specI cluster
    is_metagenomic - True or False
    geographic_location - Country, Continent, Ocean, etc.
    latitude, longitude - geographical coordinates in decimals
    general_envo_name - environment classification
    environment_material - environment material, e.g.: stool

3. **data/SPHERE_v.2022-03.levels_assessment.tsv.gz:**

Table from AMPSphere v2022-03 version. It consists in a table with the clustering 
information in three different levels (100%, 100-85%, 100-85-75%) after conversion 
using an 8-letter alphabet. Its columns are:

    AMP accession - accession in AMPSphere v2022-03 resource
    evaluation vs. representative - summary of clustering procedure
    SPHERE_fam level I - cluster accession at 1st level (100% identity)
    SPHERE_fam level II - cluster accession at 2nd level (100-85% identity)
    SPHERE_fam level III - cluster accession at 3rd level (100-85-75% identity)

4. **data/quality_families.txt.gz:**

Table generated during the quality assessment. It contains the number of clusters
passing all the tests and which percent it represents for the cluster, also
contains the information about experimental evidence of at least 1 candidate in the
cluster. Its columns are:

    family - AMP cluster accession at 3rd level (100-85-75% identity)
    experimental_evidence - True or False('*')
    quality_candidates - number of peptides passing in all quality tests ('**')
    total - total number of peptides in the cluster
    perc - percent of candidates passing in all tests

'*'  regarding to the identification of the AMP in metaproteomes or metatranscriptomes
   at least 1 candidate needed to become True

'**' simulatenously in RNAcode, AntiFam, and Terminal placement

5. **data/quality_candidates.txt.gz:**

List of candidates passing simultaneously in all quality tests (RNAcode, AntiFam, and
Terminal Placement), but without experimental evidence of translation or transcription.
This list is generated during the quality assessment and can be reobtained by using
the file `quality_assessment.tsv.gz`.

6. **data/high_quality_candidates.txt.gz:**

List of candidates passing simultaneously in all quality tests (RNAcode, AntiFam, and
Terminal Placement) with experimental evidence of translation or transcription.
This list is generated during the quality assessment and can be reobtained by using
the file `quality_assessment.tsv.gz`.

7. **data/output_clustering_significance_levelI.tsv.gz:**

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


8. **data/output_clustering_significance_levelII.tsv.gz**

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

9. **data/output_clustering_significance_levelIII.tsv.gz**

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
    
10. **data/bac120_r202.tre**

Tree of representants of species and taxa from GTDB in Newick format. It is a fixed
version from the release 202. You can find this in the following
[address](https://data.ace.uq.edu.au/public/gtdb/data/releases/release202/202.0/bac120_r202.tree).

11. **data/merged_to_donor_samples.tsv.gz**

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

12. **data/merged_to_recipient_samples_pre.tsv.gz**

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

13. **data/merged_to_recipient_samples_post.tsv.gz**

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

14. **data/dramp_v3.features.tsv.gz:**

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

15. **data/macrel_trainpos.features.tsv.gz:**

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

16. **data/ampsphere_v2022-03.features.tsv.gz:**

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

17. **data/AMPSphere_v.2022-03.annotation.tsv.gz:**

Table generated during annotation of motifs. The general regex expressions were searched
using the function __str__.contains() and results are shown by presence/absence.
Its columns are:

    id - accession in AMPSphere v2022-03 resource
    AA_rich - most frequent residue (abundance %)
    AA_absent - absent amino acids in the peptide
    sol.1 - does not start or end with a charged residue
    sol.2 - does not contain more than 1 P or G in the sequence
    sol.3 - is not composed more than 25% by a single residue
    sol.4 - no more than 45% of charged/hydrophobic residues
    sol.5 - no more than 75% of gel-prone residues
    sol.6 - net charge <= 1
    syn.1 - no forbiden motifs: 2 consecutive prolines,
            DG or DP, or N/Q at amino-terminal position
    syn.2 - charged residues <= 50% of peptide
    syn.3 - no residues that are oxidized (CMW)
    motif_match - str, contains the list of motifs per AMP ('*')
    
'*' motifs searched were available in [Ruhanen et al. (2014)](https://pubmed.ncbi.nlm.nih.gov/24478765/)
    and [Huan et al. (2021)](https://pubmed.ncbi.nlm.nih.gov/33178164/).

18. **data/freeze.v2.motusv2_5.mg3.insertcount.tsv.xz:**

Table containing the abundance of mOTUs calculated from the metagenome assemblies.
It was calculated using the mOTUs pipeline implemented in NGLess. The file is structured
with as follows:

    62,333 rows representing samples (1st column)
    14,213 columns representing mOTUs 
    
19. **data/complete_amps_associated_taxonomy.tsv.gz:**

Table obtained during metadata processing. It contains detailed information about 
AMP genes, such as their start/end codons, GC contents, GTDB annotation of their 
contigs using MMSeqs2, and etc. Its columns are:

    amp - accession in AMPSphere v2022-03 resource
    gmsc - accession in the GMSC10 resource
    sample - biosample code
    contig - contig name
    start - gene start position in the contig
    stop - gene stop position in the contig
    strand - +1 or -1 for direct or reverse complement, respectively
    fid - random unique code for the gene in the contig
    partial - ends classification (complete genes - 00)
    start_type - start codon
    rbs_motif - ribosome binding motif
    rbs_spacer - ribosome binding motif spacer
    gc_cont - content of GC in percent
    taxid - NCBI taxID
    level - taxonomy level (Root, Kingdom, ..., genus, species)
    source - microbial source name
    retained - float
    assigned - float
    agreement - float
    support - probability, support measure, float
    specI - species cluster from ProGenomes2

20. **data/bps-per-taxon.tsv.xz:**

Table containing the number of assembled base pairs per taxon
in the AMPSphere pre-analysis. Its columns are:

    taxid - NCBI taxID
    level - taxonomy level (Root, Kingdom, ..., genus, species)
    name - microbial taxon 
    nbps - number of assembled base pairs
    
21. **data/reduced_metadata.tsv.gz:**

Reduced version of metadata table. Its columns are:

    sample_accession - sample code from NCBI/ENA
    geographic_location - geographic location as country, ocean, and etc.
    latitude - decimal latitude coordinate
    longitude - decimal longitude coordinate
    general_envo_name - habitat
    environment_material - sampled material

22. **data/result_gmgc.m8.xz:**

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

23. **data/samples-min500k-assembly-prodigal-stats.tsv.gz:**

Table from the resource used to build AMPSphere v2022-03 version.
It consists of the following columns:

    sample_accession - biosample
    study - study name
    study_accession - accession of the project in NCBI/ENA
    source - source of the sample
    human_filtering - filtering off human genome reads
    inserts_raw - total number of inserts
    inserts_hq - high-quality inserts
    inserts_filtered - number of inserts after filtering by reference
    assembly_total_length - total base pairs assembled
    assembly_number - number of contigs in the assembly
    assembly_mean_length - mean contig length
    assembly_longest - longest contig length in the assembly
    assembly_shortest - shortest contig length in the assembly
    assembly_N_count - counting Ns
    assembly_Gaps - gaps in assembly
    assembly_N50 - N50 of assembly discounting Ns
    assembly_N50n - N50 of assembly counting Ns
    assembly_N70 - N70 of assembly discounting Ns
    assembly_N70n - N70 of assembly counting Ns
    assembly_N90 - N90 of assembly discounting Ns
    assembly_N90n - N90 of assembly counting Ns
    prodigal_complete - number of ORFs with start or end codons
    prodigal_no_end - number of ORFs without end codon
    prodigal_no_start - number of ORFs without start codon
    prodigal_no_start_end - number of ORFs without start or end codons
    prodigal_smorfs - number of small ORFs predicted with Prodigal
    prodigal_total_orfs - total number of ORFs
    smORFs - number of complete small ORFs

'*' Table available in this [link](https://github.com/BigDataBiology/global_data/tree/master/freeze.v2)

24. **data/metadata.tsv.gz:**

Table containing metadata description of the samples used in AMPSphere v2022-03 version.
This table can be found in the [link](https://github.com/BigDataBiology/global_data/tree/master/freeze.v2).
It consists of the following columns:

    sample_accession - biosample accession
    ena_ers_sample_id - name of the sample
    database - database available 
    access_status - status as 'public' or not
    study - study name
    study_accession - accession of study 
    publications - reference
    collection_date - date
    aliases - alternative accession
    microontology - microontology classification
    environment_biome - biome ontology
    environment_featureenvironment_material - sampled material, e.g. stool
    geographic_location - location of the sample: country, continent, ocean...
    latitude - decimal latitude
    longitude - decimal longitude
    tax_id - NCBI taxID of organism
    tax_scientific_name - taxonomy of organism
    host_common_name - host common name
    host_scientific_name - host scientific name
    host_tax_id - NCBI taxID of host
    gender - host gender
    age_years - host age
    subject_disease_status - details about disease
    antibiotic - use of antibiotics, details
    notes - commentary

25. **data/general_envo_names.tsv.gz:**

Table containing conversion tools to transform the microontology into
72 environments general habitats:

    # samples - number of samples
    # amps - number of amps detected in that habitat
    microontology - microontology term
    host_scientific_name - host species name
    host_tax_id - host NCBI taxID
    general_envo_name - general habitat name

26. **data/dramp_candidates.txt.gz**

List of AMP candidates matching to DRAMP v3 database with a significant E-value.
Based in MMSeqs2 results.

27. **data/gmgc_candidates.txt.gz**

List of AMP candidates matching to GMGC v1 database with a significant E-value.
These candidates were pre-filtered by using only AMP candidates passing the
terminal placement test - non-fragments.
Based in MMSeqs2 results.

28. **data/SmProt_candidates.txt.gz**

List of AMP candidates matching to SmProt2 database with significant E-value.
Based in MMSeqs2 results.

29. **data/starPepDB_candidates.txt.gz**

List of AMP candidates matching to starPep45k database with significant E-value.
Based in MMSeqs2 results.

30. **data/STsORFs_candidates.txt.gz**

List of AMP candidates matching to STsORFs database with significant E-value.
STsORFs are high confidence small ORFs from *Salmonella* genus.
Based in MMSeqs2 results.

31. **data/quality_assessment.tsv.gz**

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

32. **data/AMPSphere_v.2022-03.fna.xz:**

Fasta file containing the nucleotide sequences for the AMP candidates in 
AMPSphere. The arrangement of this fasta file has headers for the sequences
following the standard:

```
>GMSC10.SMORF.000_000_000_949 | AMP10.000_004
    |_ GMSC accession                |_ AMP accession
```

33. **data/bac120_taxonomy_r207.tsv.xz**

Table with the taxonomy lineages from GTDB.
It was obtained from the online resource, excluding
the first column relative to the genomes, substituting
the separatos to be tab and cleaning the duplicates.
Its columns are:

    genome - from GTDB release
    domain
    phylum
    class
    order
    family
    genus
    species

34. **data/adjust_significant_function.csv.xz**

Enrichment analysis of functions among those ortholog groups identified by mapping
AMPs against GMGC. This analysis was done using all proteins in GMGC not only the
unigenes, which were deduplicated prior counting. Its columns are as follow:

    eggnog_OG - EggNOG ortholog group
    count_AMP - number of AMPs matching to an OG
    count_GMGC - number of proteins from GMGC matching an OG
    total_AMP - total AMPs matched against GMGC
    total_GMGC - total proteins in GMGC
    p_value - p_value calculated from the hypergeometric space
    amp_fraction - proportion of total AMPs matching an OG
    GMGC_fraction - proportion of large proteins matching an OG
    times - fold enrichment (AMP fraction / GMGC fraction)
    p_adjust - p-value adjusted using Bonferroni
    sig - if p<0.05
    label - decoy

35. **data/amp_COG.tsv.xz**

Table contaning the COG classes for the mapped AMPs. It does contain about 6k entries
because most of the OGs from GMGC were not in the metadata table relative to the 
EggNOG mapping. Its columns are:

    cog_class - class from [COG database](https://ecoliwiki.org/colipedia/index.php/Clusters_of_Orthologous_Groups_(COGs))
    count_AMP - mapped AMPs
    fraction - proportion of mapped AMPs


