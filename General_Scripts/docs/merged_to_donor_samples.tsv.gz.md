**merged_to_donor_samples.tsv**

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

**MD5 SUM:**	df6d568bc12b57d8245f0bc52ef00ebd

**Size (MBytes):**	0.04468727111816406

**Content sample (first 5 items):**

species	fmt.id	fmt.type	timepoint.fmt_x	sample.donor	sample.pre	sample.post	outcome	sp_retained	reject_consp	coexist_consp.rec	coexist_consp.donor	engraft_consp	reject_novel	engraft_novel	sp_lost	influx_consp	influx_novel	influx_total	genome	completeness	contamination	subject_type	timepoint.fmt_y	clinical_response	ORFs	smORFs	#_amps	sequences
ANI_AL_95_00140	FAME.06b	FMT	7.0	EBS_DON_014-11-0-0	EBS_014-11-0-0	EBS_014-11-7-0	rejection novel	0.0	0.0	0.0	0.0	0.0	1.0	0.0	0.0	0.0	0.0	0.0	FAME.per_sample.EBS_DON_014-11-0-0_16s007909-1-2_lane2.screened.adapter.screened.hg19.0026	92.62	0.0	donor	0		2375.0	217.0	1.0	SMKTLSVFMFSHVHLKRLIYEYFKI
ANI_AL_95_00140	FAME.06b	FMT	14.0	EBS_DON_014-11-0-0	EBS_014-11-0-0	EBS_014-11-14-0	novel unclear	0.0	0.0	0.0	0.0	0.0	0.0	0.11511627906976744	0.0	0.0	0.8848837209302326	0.8848837209302326	FAME.per_sample.EBS_DON_014-11-0-0_16s007909-1-2_lane2.screened.adapter.screened.hg19.0026	92.62	0.0	donor	0		2375.0	217.0	1.0	SMKTLSVFMFSHVHLKRLIYEYFKI
ANI_AL_95_00140	FAME.10	FMT	7.0	EBS_DON_015-11-0-0	EBS_015-11-0-0	EBS_015-11-7-0	rejection conspecific	0.0	0.7733522996680892	0.0	0.08036984352773827	0.0	0.0	0.0	0.0	0.14627785680417257	0.0	0.14627785680417257	FAME.per_sample.EBS_DON_015-11-0-0_16s007917-1-1_lane8.screened.adapter.screened.hg19.0049	98.66	0.0	donor	0		2387.0	216.0	1.0	SMKTLSVFMFSHVHLKRLIYEYFKI
ANI_AL_95_00140	FAME.10	FMT	28.0	EBS_DON_015-11-0-0	EBS_015-11-0-0	EBS_015-11-28-0	species lost	0.0	0.0	0.0	0.0	0.0	0.0	0.0	1.0	0.0	0.0	0.0	FAME.per_sample.EBS_DON_015-11-0-0_16s007917-1-1_lane8.screened.adapter.screened.hg19.0049	98.66	0.0	donor	0		2387.0	216.0	1.0	SMKTLSVFMFSHVHLKRLIYEYFKI
[...]
