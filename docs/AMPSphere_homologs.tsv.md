**AMPSphere_homologs**

**Description:**	gives the homology tests against DRAMP, SmProt and StarPep45k
  
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

**MD5 SUM:**	3db754f116cfb54d95fff249517b33e7

**Size (MBytes):**	100.63867664337158

**Content sample (first 5 items):**

query	target	evalue	gapopen	pident	nident	qstart	qend	qlen	tstart	tend	tlen	alnlen	raw	bits	cigar	qseq	tseq	qheader	theaderqaln	taln	qframe	tframe	mismatch	qcov	tcov
AMP10.292_996	DRAMP14850_PATENTED	8.169e-11	0	81.4	22	127	30	58	84	86	27	113	49	27M	ITPRGEKKAYVKLKPEYKASELAVKLGVIV	MDPYKVIIRPVVTEKAISLIEKENKLTFIVDRRATKTDIKKAIEEIFNVKVEKVNTLITPKGEKKAYVKLKPEYSASEIAARLGLF	AMP10.292_996 | SPHERE-III.016_390	DRAMP14850_PATENTED	ITPRGEKKAYVKLKPEYKASELAVKLG	ITPKGEKKAYVKLKPEYSASEIAARLG	50.9	0.314
AMP10.292_996	DRAMP14997_PATENTED	4.904e-08	0	67.8	19	128	30	58	85	86	28	95	42	28M	ITPRGEKKAYVKLKPEYKASELAVKLGVIV	MDAFDVIKAPVVTEKTVRMIEEENKLVFYVDRRATKQDIKRAMKELFDVEVEKVNTLITPKGEKKAYVKLKEGYDASKIAASLGIY	AMP10.292_996 | SPHERE-III.016_390	DRAMP14997_PATENTED	ITPRGEKKAYVKLKPEYKASELAVKLGV	ITPKGEKKAYVKLKEGYDASKIAASLGI	90.933	0.326
AMP10.617_764	DRAMP15338_PATENTED	3.641e-07	0	52.9	18	336	37	54	87	88	34	90	40	34M	ALAEPELMQAAQRNIIHKNNASRKVSRLAHQIAKLAK	MANTTSAKKATRKIARRSAVNKARRSRIRSFVRKVEEAIASGDQALAAAALKAAQPELMRPATKGVMHSNTASRKVSRLAQRVKSLSA	AMP10.617_764 | SPHERE-III.001_940	DRAMP15338_PATENTED	AEPELMQAAQRNIIHKNNASRKVSRLAHQIAKLA	AQPELMRPATKGVMHSNTASRKVSRLAQRVKSLS			16	0.919	0.386
AMP10.207_399	DRAMP15280_PATENTED	2.582e-07	0	61.2	19	131	31	59	89	89	31	90	40	31M	LVGRRRRFLNYLQKKNLEGYRSLIRELGLRR	MALTAEQKKEILRSYGLHETDTGSPEAQIALLTKRIADLTEHLKVHKHDHHSRRGLLLLVGRRRRLIKYISQIDVERYRSLIERLGLRR	AMP10.207_399 | SPHERE-III.000_004	DRAMP15280_PATENTED	LVGRRRRFLNYLQKKNLEGYRSLIRELGLRR	LVGRRRRLIKYISQIDVERYRSLIERLGLRR			12	1.0	0.348
[...]
