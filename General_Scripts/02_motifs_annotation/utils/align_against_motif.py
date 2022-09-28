def align_against_motif(x: str) -> str:
    '''
    Align protein X against motif
       Regex rule: [KR]{0,2}[KR].{0,2}[KR].{2,4}[ILVM].[ILVF]

    Input: str
    Output: str
    '''
    x = list(x)
    while (len(x) < 9) or (x[-8] not in ['K', 'R']):
        x.insert(-3, '-')
    while len(x[:-8]) < 5:
        if len(x[:-8]) == 1:
            x[:-8] = ['-', '-', x[:-8][0], '-', '-']
        else:
            while len(x[:-8]) < 3: 
                x.insert(-8, '-')
            while x[:-8][-3] not in ['K', 'R']:
                x.insert(-8, '-')
            while len(x) < 13: x.insert(0, '-')
    return ''.join(x)
