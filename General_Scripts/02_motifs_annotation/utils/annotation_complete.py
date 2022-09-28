import re
import numpy as np
from Bio import SeqIO


def entropy(labels):
    prob_dict = {x:labels.count(x)/len(labels) for x in labels}
    probs = np.array(list(prob_dict.values()))
    return -probs.dot(np.log2(probs))


def rich_poor(labels):
    aa = ['A','C','D','E','F',
          'G','H','I','K','L',
          'M','N','P','Q','R',
          'S','T','V','W','Y']
    prob_dict = {x:labels.count(x)*100/len(labels) for x in labels}
    key, value = max(prob_dict.items(), key=lambda x:x[1])
    absent_aa = [ x for x in aa if x not in prob_dict]    
    return f'{key} ({value:.1f}%)', absent_aa


def gene_probability(labels):
    aa = {'A':4/64,'C':2/64,'D':2/64,'E':2/64,'F':2/64,
          'G':4/64,'H':2/64,'I':3/64,'K':2/64,'L':6/64,
          'M':1/64,'N':2/64,'P':4/64,'Q':2/64,'R':6/64,
          'S':6/64,'T':4/64,'V':4/64,'W':1/64,'Y':2/64,
          '*':3/64}
    prob = np.prod([aa[x] for x in labels])
    return prob


def solubility_rule1(sequence):
    charged_residues = 'HRKDE'
    if sequence[0] in charged_residues:
        return False
    elif sequence[-1] in charged_residues:
        return False
    else:
        return True
   

def solubility_rule2(sequence):
    P = sequence.count('P')
    G = sequence.count('G')
    if P+G > 1:
        return False
    else:
        return True


def solubility_rule3(sequence):
    from collections import Counter
    d = dict(Counter(sequence))
    suma = sum(d.values())
    for k, v in d.items():
        if v/suma > 0.25:
            return False
        else:
            return True


def solubility_rule4(sequence):
    charged = 'HRKDE'
    hydrophob = 'VILMFWC'
    res = 0
    for cha in sequence:
        if cha in charged or cha in hydrophob:
            res += 1
    if res/len(sequence) > 0.45:
        return False
    else:
        return True


def solubility_rule5(sequence):
    gel_prone = 'DEHKNQRSTY'
    res = 0
    for cha in sequence:
        if cha in gel_prone:
            res += 1
    if res/len(sequence) > 0.75:
        return False
    else:
        return True
        

def solubility_rule6(sequence):
    pos = 'KRH'
    neg = 'DE'
    net = 0
    for i in sequence:
        if i in pos: 
            net += 1
        elif i in neg:
            net -= 1
    if net > 1:
        return False
    else:
        return True


def synthesis_rule1(sequence):
    import re
    forbidden_motifs = {'2-prolines': r'[P]{3,}', 'DG-DP': r'D[GP]', 'N-Q-Nterminal': r'^[NQ]'}
    for motif in forbidden_motifs:
        if re.search(forbidden_motifs[motif], sequence):
            return False
    return True


def synthesis_rule2(sequence):
    charged_residues = 'HRKDE'
    counter_charged = 0
    cc = []
    for residue in sequence:
         counter_charged += 1
         if residue in charged_residues:
            counter_charged = 0
         if counter_charged >= 5:
            cc.append(1)
         else:
            cc.append(0)
    t = cc.count(0) / len(cc)
    if t > 0.5:
        return False
    else:
        return True


def synthesis_rule3(sequence):
    aa_oxidation = 'MCW'
    for cha in sequence:
        if cha in aa_oxidation:
            return False
    return True


def load_motif():
    import pandas as pd
    data = pd.read_table('data/db_motif.tsv',
                         sep='\t',
                         header='infer')
    dbmotif = data[['motif', 'name']]
    dbmotif = dbmotif.groupby('name')
    dbmotif = dbmotif.agg(lambda x: list(x))
    return dbmotif.to_dict()['motif']    
    
    
def anno_motif(sequence, dbmotif):
    motif_match = []
    for name, mlist in dbmotif.items():
        for motif in mlist:
            if re.search(motif, sequence):
                motif_match.append(name)
    return set(motif_match)


def pipe_anno(infile, ofile):
    dbmotif = load_motif()
    with open(ofile, 'w') as fout:
        fout.write('id\tAA_rich\tAA_absent\tH\tlog2(gene_prob)\tsol.1\tsol.2\tsol.3\tsol.4\tsol.5\tsol.6\tsyn.1\tsyn.2\tsyn.3\tmotif_match\n')
        for record in SeqIO.parse(infile, 'fasta'):
            record.seq = str(record.seq)
            H = entropy(record.seq)
            rich, poor = rich_poor(record.seq)
            poor = '|'.join(poor)
            prob = gene_probability(record.seq)
            sol1 = solubility_rule1(record.seq)
            sol2 = solubility_rule2(record.seq)
            sol3 = solubility_rule3(record.seq)
            sol4 = solubility_rule4(record.seq)
            sol5 = solubility_rule5(record.seq)
            sol6 = solubility_rule6(record.seq)
            syn1 = synthesis_rule1(record.seq)
            syn2 = synthesis_rule2(record.seq)
            syn3 = synthesis_rule3(record.seq)
            motif_match = anno_motif(record.seq, dbmotif)
            motif_match = '|'.join(motif_match)
            fout.write(f'{record.id}\t{rich}\t{poor}\t{H}\t{np.log2(prob)}\t{sol1}\t{sol2}\t{sol3}\t{sol4}\t{sol5}\t{sol6}\t{syn1}\t{syn2}\t{syn3}\t{motif_match}\n')
            print(f'{record.id}')
            
