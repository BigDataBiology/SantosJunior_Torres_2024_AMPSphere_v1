**AMPSphere_v.2022-03.annotation.tsv**

**Description:**	Table generated during annotation of motifs. The general regex expressions were searched
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

**MD5 SUM:**	97d5e310a8707ec54092735c6738ca18

**Size (MBytes):**	9.27689266204834

**Content sample (first 5 items):**

id	AA_rich	AA_absent	sol.1	sol.2	sol.3	sol.4	sol.5	sol.6	syn.1	syn.2	syn.3	motif_match
AMP10.000_000	K (25.9%)	C|D|H|P|Q|R|T|Y	False	False	False	False	True	False	True	False	False	gala|adherence|capsule|lipid binding motif
AMP10.000_001	R (14.0%)	C|H|P|W	True	False	True	False	True	False	True	False	False	gala|gly zipper|autoproteolytic cleavage motif
AMP10.000_002	K (14.8%)	C|D|H|Q|T|W	False	False	True	False	True	False	True	False	False	
AMP10.000_003	V (20.6%)	C|F|H|N|P|S|T|W|Y	False	False	True	False	True	False	True	False	False	dock_MAPK|gly zipper
[...]
