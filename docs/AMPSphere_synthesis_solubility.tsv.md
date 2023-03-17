**AMPSphere_synthesis_solubility**

**Description:**	gives information about the criteria the AMPs
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

**MD5 SUM:**	292859e422355f8ac6e5997484870b2d

**Size (MBytes):**	52.764655113220215

**Content sample (first 5 items):**

amp	sol_rule_1	sol_rule_2	sol_rule_3	sol_rule_4	sol_rule_5	sol_rule_6	syn_rule_1	syn_rule_2	syn_rule_3
AMP10.000_000	False	False	False	False	True	False	True	False	False
AMP10.000_001	True	False	True	False	True	False	True	False	False
AMP10.000_002	False	False	True	False	True	False	True	False	False
AMP10.000_003	False	False	True	False	True	False	True	False	False
[...]
