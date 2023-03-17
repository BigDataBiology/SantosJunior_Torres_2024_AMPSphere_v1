**merged_to_recipient_samples_pre.tsv**

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



**MD5 SUM:**	f8000671120e7e22f6cd60b97e261a5e

**Size (MBytes):**	0.05611991882324219

**Content sample (first 5 items):**

species	fmt.id	fmt.type	timepoint.fmt_x	sample.donor	sample.pre	sample.post	outcome	sp_retained	reject_consp	coexist_consp.rec	coexist_consp.donor	engraft_consp	reject_novel	engraft_novel	sp_lost	influx_consp	influx_novel	influx_total	genome	completeness	contamination	subject_type	timepoint.fmt_y	clinical_response	ORFs	smORFs	#_amps	sequences
ANI_AL_95_00294	FATLOSE.11	FMT	21.0	FATLOSE_DON_002-22-0-0	FATLOSE002-22-0-0	FATLOSE002-22-21-0	coexistence	0.0	0.0	0.1201019664967225	0.7833272638990046	0.0	0.0	0.0	0.0	0.09657076960427291	0.0	0.09657076960427291	FATLOSE.per_sample.FATLOSE002-22-0-0.0002	94.0	0.0	recipient	0		2071.0	254.0	2.0	KPPSFLNWSINLYTCNTESCQTLGQLNRMKGSKKMYQ,FTVFILPFRFVHSVKKQDSALKTRRALFHFFVYVL
ANI_AL_95_00294	FATLOSE.11	FMT	84.0	FATLOSE_DON_002-22-0-0	FATLOSE002-22-0-0	FATLOSE002-22-84-0	coexistence	0.0	0.0	0.24331488693038192	0.6400765184122429	0.0	0.0	0.0	0.0	0.11660859465737516	0.0	0.11660859465737516	FATLOSE.per_sample.FATLOSE002-22-0-0.0002	94.0	0.0	recipient	0		2071.0	254.0	2.0	KPPSFLNWSINLYTCNTESCQTLGQLNRMKGSKKMYQ,FTVFILPFRFVHSVKKQDSALKTRRALFHFFVYVL
ANI_AL_95_00294	FATLOSE.11	FMT	126.0	FATLOSE_DON_002-22-0-0	FATLOSE002-22-0-0	FATLOSE002-22-126-0	coexistence	0.0	0.0	0.20619127988748243	0.6275105485232068	0.0	0.0	0.0	0.0	0.16629817158931082	0.0	0.16629817158931082	FATLOSE.per_sample.FATLOSE002-22-0-0.0002	94.0	0.0	recipient	0		2071.0	254.0	2.0	KPPSFLNWSINLYTCNTESCQTLGQLNRMKGSKKMYQ,FTVFILPFRFVHSVKKQDSALKTRRALFHFFVYVL
ANI_AL_95_00294	FAME.06b	FMT	7.0	EBS_DON_014-11-0-0	EBS_014-11-0-0	EBS_014-11-7-0	species lost	0.0	0.0	0.0	0.0	0.0	0.0	0.0	1.0	0.0	0.0	0.0	FAME.per_sample.EBS_014-11-0-0_17s000833-1-1_lane5.screened.adapter.screened.hg19.0029	93.96	0.0	recipient	0	non-responder	1951.0	230.0	0.0	
[...]
