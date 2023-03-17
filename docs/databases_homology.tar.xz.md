**databases_homology.tar**

**Description:**	Databases used for homologs detection. It is structured as a tar.xzipped folder.
                        The folder contains different databases and the pre-computed results for the 
                        search of AMPSphere against the GMGC v1 dataset.

Databases included as fasta files with peptide sequences:

  (i) [DRAMP v.3](dramp.cpu-bioinfor.org/)
  (ii) [starPepDB 45k](http://mobiosd-hub.com/starpep)
  (iii) [smProt](bigdata.ibp.ac.cn/SmProt)
  (iv) [STsORFs](https://academic.oup.com/microlife/article/1/1/uqaa002/5928550)
  
Databases included as pre-computed resources:

  (v) true_pep_2022_vs_progenomesgmgc.tsv.xz - output obtained from the MMSeqs2 easy-search algorithm
                                               using the peptide sequences from AMPSphere as query and
                                               the GMGC v1 database (including ProGenomes sequences)
                                               as the target. The output includes the following output
                                               format: query, target, evalue, gapopen, pident, nident,
                                               qstart, qend, qlen, tstart, tend, tlen, alnlen,
                                               raw, bits, cigar, qseq, tseq, qheader, theader, qaln,
                                               taln, qframe, tframe, mismatch, qcov, tcov
  
**MD5 SUM:**	1491346cb6f299296efce1dfc7507393

**Size (MBytes):**	1110.340877532959

**Content sample (first 5 items):**

databases_homology
databases_homology/DRAMP.fa
databases_homology/converter
databases_homology/converter/AMPSphere_2021_to_2022.py
databases_homology/converter/data
databases_homology/converter/data/AMPSphere_v.2022-03.faa.gz
databases_homology/converter/data/true_pep_vs_progenomesgmgc.tsv.xz
databases_homology/converter/data/AMPSphere_v.2021-03.faa.gz
databases_homology/converter/README.md
databases_homology/converter/analysis
databases_homology/converter/analysis/true_pep_2022_vs_progenomesgmgc.tsv.xz
databases_homology/converter/analysis/converter_AMPSphere_2021_2022.py
databases_homology/starPepDB.fasta
databases_homology/true_pep_2022_vs_progenomesgmgc.tsv.xz
databases_homology/all_SmProt.fa.gz
databases_homology/STsORFs.faa
