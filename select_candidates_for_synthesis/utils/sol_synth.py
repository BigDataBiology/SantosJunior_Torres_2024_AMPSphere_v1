class solubility:

    sequence: str
    
    def __init__(self, sequence):
        self.sequence = sequence
    

    def rule1(self):
        charged_residues = 'HRKDE'
        if self.sequence[0] in charged_residues:
            return False
        elif self.sequence[-1] in charged_residues:
            return False
        else:
            return True
   

    def rule2(self):
        P = self.sequence.count('P')
        G = self.sequence.count('G')
        if P+G > 1:
            return False
        else:
            return True


    def rule3(self):
        from collections import Counter
        d = dict(Counter(self.sequence))
        suma = sum(d.values())
        for k, v in d.items():
            if v/suma > 0.25:
                return False
            else:
                return True


    def rule4(self):
        charged = 'HRKDE'
        hydrophob = 'VILMFWC'
        res = 0
        for cha in self.sequence:
            if cha in charged or cha in hydrophob:
                res += 1
        if res/len(self.sequence) > 0.45:
            return False
        else:
            return True


    def rule5(self):
        gel_prone = 'DEHKNQRSTY'
        res = 0
        for cha in self.sequence:
            if cha in gel_prone:
                res += 1
        if res/len(self.sequence) > 0.75:
            return False
        else:
            return True
            

    def rule6(self):
        pos = 'KRH'
        neg = 'DE'
        net = 0
        for i in self.sequence:
            if i in pos: 
                net += 1
            elif i in neg:
                net -= 1
        if net > 1:
            return False
        else:
            return True

            
    def rule_judgement(self):
        r1 = self.rule1()
        r2 = self.rule2()
        r3 = self.rule3()
        r4 = self.rule4()
        r5 = self.rule5()
        r6 = self.rule6()
        return [r1, r2, r3,
                r4, r5, r6]
                
                
class synthesis:

    sequence: str
    
    def __init__(self, sequence):
        self.sequence = sequence
        
    
    def rule1(self):
        import re
        forbidden_motifs = {'2-prolines': r'[P]{3,}', 'DG-DP': r'D[GP]', 'N-Q-Nterminal': r'^[NQ]'}
        for motif in forbidden_motifs:
            if re.search(forbidden_motifs[motif], self.sequence):
                return False
        return True

    
    def rule2(self):
        charged_residues = 'HRKDE'
        counter_charged = 0
        cc = []
        for residue in self.sequence:
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


    def rule3(self):
        aa_oxidation = 'MCW'
        for cha in self.sequence:
            if cha in aa_oxidation:
                return False
        return True
    
   
    def rule_judgement(self):
        r1 = self.rule1()
        r2 = self.rule2()
        r3 = self.rule3()
        return [r1, r2, r3]        
        
