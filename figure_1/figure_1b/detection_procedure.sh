format="query,target,evalue,gapopen,pident,nident,qstart,qend,qlen,tstart,tend,tlen,alnlen,raw,bits,cigar,qseq,tseq,qheader,theader,qaln,taln,qframe,tframe,mismatch,qcov,tcov"

# AMPSphere_v.2021-03.faa.gz is available in the Zenodo repository
# databases used in the mmseqs searches are in ubuntu@aws.big-data-biology.org:/share/work/Celio/files_for_figures/databases_homology/

# search for homologs
mkdir -p DRAMP SmProt STsORFs starpep GMGC10_pGenomes_db
mmseqs easy-search AMPSphere_v.2021-03.faa.gz DRAMP.fa DRAMP/result_dramp.m8 tmp --format-output $format
mmseqs easy-search AMPSphere_v.2021-03.faa.gz all_SmProt.fa.gz SmProt/result_SmProt.m8 tmp --format-output $format
mmseqs easy-search AMPSphere_v.2021-03.faa.gz starPepDB.fasta starpep/result_starPepDB.m8 tmp --format-output $format
mmseqs easy-search AMPSphere_v.2021-03.faa.gz STsORFs.faa STsORFs/result_STsORFs.m8 tmp --format-output $format

# special search in a large database (GMGCv1)
python3 splitGMGC10pgen.py

for chunk in $(ls GMGC10_pGenomes_db/*.fasta):
do
	mmseqs easy-search AMPSphere_v.2021-03.faa.gz $chunk ${chunk/.fasta/.m8} tmp --disk-space-limit 5G \
                                      --db-load-mode 3 --format-output $format --threads 3
done

cat GMGC10_pGenomes_db/*m8 | awk '$3 <= 1e-5' | xz > truepep_gmgc_progenomes.m8.xz

# retrieving statistically significant matches
awk '$3 <= 1e-5 {print $1}' DRAMP/result_dramp.m8 | sort | uniq > dramp_candidates.txt
awk '$3 <= 1e-5 {print $1}' SmProt/results/most_recent_all/result_SmProt.m8 | sort | uniq > smprots_candidates.txt
awk '$3 <= 1e-5 {print $1}' starpep/result_starPepDB.m8 | sort | uniq > star_candidates.txt
awk '$3 <= 1e-5 {print $1}' STsORFs/result_STsORFs.m8 | sort | uniq > STsORFs_candidates.txt
xzcat GMGC/truepep_gmgc_progenomes.m8.xz | awk '$3 <= 1e-5 {print $1}' | sort | uniq > gmgc_candidates.txt

## check number of homologs by database:
wc -l dramp_candidates.txt  ## detected in DRAMP3 (6,339)
wc -l smprots_candidates.txt  ## detected in SmProt2 (15,297)
wc -l star_candidates.txt  ## detected in starPep45k (6,263)
wc -l STsORFs_candidates.txt  ## detected in STsORFs (11)
wc -l gmgc_candidates.txt  ## detected in GMGCv1 (61,020)

## make all homologs in one file
cat *_candidates.txt | sort | uniq > annotated_candidates.txt

## getting quality candidates
# files: high_quality_candidates.txt and quality_candidates.txt
# available in ubuntu@aws.big-data-biology.org:/share/work/Celio/files_for_figures/quality_control
cat high_quality_candidates.txt quality_candidates.txt | sort | uniq > qc_controlled_amps.txt

# checking overlaps
comm -1 -2 annotated_candidates.txt qc_controlled_amps.txt  | wc -l  ## detected in both (10427)
comm -3 -2 annotated_candidates.txt qc_controlled_amps.txt  | wc -l  ## detected in annotated only (63347)
comm -3 -1 annotated_candidates.txt qc_controlled_amps.txt  | wc -l  ## detected in qc_controlled only (19769)
