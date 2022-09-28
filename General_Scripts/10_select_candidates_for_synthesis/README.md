# AMPSphere v2022_03

The following files structure the AMPSphere resource:

I) AMPSphere_features.tsv: list of features predicted from the primary sequence of peptides.

  Columns:
  amp - access code in AMPSphere
  charge - charge of the peptide at pH 7.0 (EMBOSS scale)
  pI - isoelectric point
  aindex - aliphatic index
  instaindex - instability index 
  boman - boman index (measures propensity to bind membranes)
  hydrophobicity - hydrophobicity by the KyteDoolitle scale
  hmoment - hydrophobic moment (measures amphiphilicity at 90 degrees of the plane 
                                established in an geometrical helix, it gives the
                                measurement of how structured charges and hydrophobic
                                residues are in the protein)

II) AMPSphere_fmt_involved.tsv: list of AMPs found in the fecal matter transplant story.

III) AMPSphere_homologs.tsv: gives the homology tests against DRAMP, SmProt and StarPep45k
  
   Peptides were searched using MMSeqs2.
   The output is a blast-like result table.

   Columns:
   query - AMPSphere peptide
   target - database matched protein
   evalue - expected value
   gapopen - gaps
   pident - percent of identity
   nident - identity counts
   qstart - start position of alignment on query
   qend - end position of alignment on query
   qlen - length of query
   tstart - start position of alignment on target
   tend - end position of alignment on target
   tlen - length of target
   alnlen - alignment length
   raw - -
   bits - bit score
   cigar - CIGAR code
   qseq - query sequence
   tseq - target sequence
   qheader - query header
   theader - target header
   qaln - query aligned segment
   taln - target aligned segment
   qframe - query frame (not used)
   tframe - target frame (not used)
   mismatch - number of mismatches
   qcov - query coverage
   tcov - target coverage

IV) AMPSphere_quality_control.tsv:  AMPs were submitted to different
                                    tests to ensure their quality.
   More info, please visit:
   [AMPSphere quality link](https://ampsphere.big-data-biology.org/quality_tests)

   Columns: 
   amp - AMP access code in AMPSphere
   Antifam - checks spurious gene predictions
   RNAcode - checks gene signatures (Passed, Failed, Not tested*)
   metaproteomes - checks for peptides in metaproteomic studies from PRIDE
   metatranscriptomes - checks for peptides mapping to at least two reads from metatranscriptomes
   terminal placement - checks if it represents a gene fragment

   * refer to those peptides that due to experimental conditions
     could not be tested


V) AMPSphere_selected_candidates_coprediction.tsv:

   Candidates passing in at least 60% of the solubility/synthesis criteria
   along with all quality criteria (about 5.6k) were assessed. These peptides
   undergone a test of AMP prediction in other 6 different systems:
   AI4AMP, amPEPpy, ampir, AMPScanner v2, APIN, AMPlify.

   These results were merged in a table with the following columns:
   amp - AMPSphere access code
   AI4AMP, AmPEPpy, ampir, AMPlify, AMPScannerv2, and APIN - prediction (AMP/non-AMP)
   pcts - percent of co-prediction as AMP
   prediction - if pcts > 0.5 then it was predicted as AMP


VI) AMPSphere_sequences.tsv: gives the sequence information for all AMPs in the 
                             AMPSphere.

    Columns:
    amp, family - access codes for the AMP and its family in AMPSphere resource
    length_residues - length of sequence in residues
    sequence - sequence of amino acids



VII) AMPSphere_synthesis_solubility.tsv: gives information about the criteria the AMPs
                                         passed or failed for solubility and synthesis.

    Solubility rules are based on recommendations of [SolyPep](http://bioserv.rpbs.univ-paris-diderot.fr/services/SolyPep/).

    Synthesis rules are based on empirical observations from the [PepFun](https://github.com/rochoa85/PepFun) group.

    Columns:
    amp - access code in the AMPSphere resource
    sol_rule_1 - not start or end with a charged residue
    sol_rule_2 - not contain more than 1 P or G in the sequence
    sol_rule_3 - not be composed more than 25% by a single residue
    sol_rule_4 - not have more than 45% of charged/hydrophobic residues
    sol_rule_5 - not have more than 75% of gel-prone residues
    sol_rule_6 - not have net charge above 1
    syn_rule_1 - not have prohibited motifs: 2 consecutive prolines, DG or DP, or N/Q at amino-terminal position
    syn_rule_2 - charged residues should not constitute more than 50% of the residue
    syn_rule_3 - no residues that are oxidized (CMW)


VIII) AMPSphere_taxonomy_environment.tsv: gives detailed information about the origin of 
                                          the different AMPs 

    Columns:
    amp - access code in the AMPSphere resource
    copies - number of samples it was detected
    source - microbial sources of a given amp (list format)
    environment - list of habitats in which it was detected
    sample - access codes of samples in which it was detected


IX) selected_candidates.tsv: table containing the selected peptides list with 364 sequences
                             For a peptide to be listed, it was:
                                  - co-predicted as AMPs by all tested methods,
                                  - passed in at least 60% of solubility and
                                    synthesis criteria,
                                  - 40 residues or shorter 

    Columns:

    amp, family - access codes in the AMPSphere resource
    copies - number of samples it was detected
    fscore - geometric mean telling the proportion of the solubility/synthesis criteria it passed
    sequence - amino acids sequence
    feature - class it belongs ('high copy', 'low copy', 'homolog', 'fmt-related')
    source - microbial sources of a given amp (list format)
    environment - list of habitats in which it was spotted

X) selected_peptides.fasta: fasta file containing suggested peptides for testing.
                            Among the set of 50 peptides selected, we included three different classes: 15 homologs matching
                            in some extension to peptides from databases such as DRAMP and smProt; 1 peptide found among
                            those in the FMT story, and the 34 most frequent peptides in our set, as "high copies".

   Format of header: >AMP10.020_784 | SPHERE-III.017_034 | 26 | high_copy
                         |_ access         |                |       |
                                           |_ family        |       |
                                                            |_ copy number
                                                                    |_ feature class
