**merged_to_recipient_samples_post.tsv**

**Description:**	Merged information of preprocessed MAGs from Schmidt et al. (2022).
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


**MD5 SUM:**	4d63c645bae3c06c963c6f2df7119388

**Size (MBytes):**	0.07717323303222656

**Content sample (first 5 items):**

species	fmt.id	fmt.type	timepoint.fmt_x	sample.donor	sample.pre	sample.post	outcome	sp_retained	reject_consp	coexist_consp.rec	coexist_consp.donor	engraft_consp	reject_novel	engraft_novel	sp_lost	influx_consp	influx_novel	influx_total	genome	completeness	contamination	subject_type	timepoint.fmt_y	clinical_response	ORFs	smORFs	#_amps	sequences
ANI_AL_95_00294	FATLOSE.17	FMT	21.0	FATLOSE_DON_026-22-0-0	FATLOSE026-22-0-0	FATLOSE026-22-21-0	coexistence	0.0	0.0	0.6958071156684534	0.13193311582381728	0.0	0.0	0.0	0.0	0.17225976850772928	0.0	0.17225976850772928	FATLOSE.per_sample.FATLOSE026-22-21-0.0070	92.09	0.0	recipient	21		1847.0	218.0	0.0	
ANI_AL_95_00333	FATLOSE.15	FMT	21.0	FATLOSE_DON_022-22-0-0	FATLOSE022-22-0-0	FATLOSE022-22-21-0	coexistence	0.0	0.0	0.6144496039333516	0.3361103523627424	0.0	0.0	0.0	0.0	0.049440043703906034	0.0	0.049440043703906034	FATLOSE.per_sample.FATLOSE022-22-21-0.0033	92.62	0.0	recipient	21		1843.0	206.0	1.0	ILGAGIYYLAKEKTDREARKIYTITTVAGALILAGAILKLVLAR
ANI_AL_95_02951	FATLOSE.22	FMT	84.0	FATLOSE_DON_048-22-0-0	FATLOSE048-22-0-0	FATLOSE048-22-84-0	novel unclear	0.0	0.0	0.0	0.0	0.0	0.0	0.14431486880466474	0.0	0.0	0.8556851311953353	0.8556851311953353	FATLOSE.per_sample.FATLOSE048-22-84-0.0036	93.51	1.57	recipient	84		2433.0	244.0	2.0	TVGSRNMMNMFAMCMGSMCMCRCANYHKVLSR,QVPGAVRFAAVRGFLQRIQNFSAIKAKKRHSLTVRW
ANI_AL_95_03045	FATLOSE.12	FMT	21.0	FATLOSE_DON_005-22-0-0	FATLOSE005-22-0-0	FATLOSE005-22-21-0	coexistence	0.0	0.0	0.779990121627462	0.11798481200222263	0.0	0.0	0.0	0.0	0.10202506637031541	0.0	0.10202506637031541	FATLOSE.per_sample.FATLOSE005-22-21-0.0094	91.69	1.34	recipient	21		1674.0	163.0	1.0	ITAVMVTVGYKVVSAVEWLMDRLLHLK
[...]
