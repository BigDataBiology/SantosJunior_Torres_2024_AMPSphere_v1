def variation(transl):
    from itertools import product
    variants = []
    for idx, l in enumerate(transl):
        if idx == 0: variants = l
        elif idx <= len(transl):
            variants = [''.join(x) for x in product(variants, l)]
    return variants
    

def hamming_distance(str1, str2):
    assert len(str1) == len(str2)
    return sum(chr1 != chr2 for chr1, chr2 in zip(str1, str2))


def mass_comp(seqlist):
    import random
    from itertools import product
    size = len(seqlist)**2
    size = int(0.1*size)
    if size > 10_000: size = 10_000
    if size < 1000: size = 1000
    identities = []
    L = len(seqlist[0])
    for _ in range(size):
        s1, s2 = 'a', 'a'
        while s1 == s2:
            s1, s2 = random.choices(seqlist, k=2)
        f = hamming_distance(s1, s2)/L
        f = 1 - f
        identities.append(f)
    return identities


def simseq(n):
    import random
    from numpy import mean, std
    minalph = {'L': ['L', 'V', 'M', 'I', 'C'],
               'A': ['A', 'G'],
               'S': ['S', 'T'],
               'P': ['P'],
               'F': ['F', 'W', 'Y'],
               'E': ['E', 'D', 'Q', 'N'],
               'K': ['K', 'R'],
               'H': ['H']}
    s = random.choices(list(minalph.keys()), weights=[5, 2, 2, 1, 3, 4, 2, 1], k=n)
    s = ''.join(s)
    transl = [minalph[i] for i in s]
    variants = variation(transl)
    variants = mass_comp(variants)
    return [s, len(transl), max(variants), min(variants), mean(variants), std(variants)]


def simulation(nseqs, length):
    import pandas as pd
    test = []
    for i in range(nseqs):
        test.append(simseq(length))
    test = pd.DataFrame(test, columns=['sequence', 'variants', 'max_id', 'min_id', 'avg_id', 'std_id'])
    return test
    

def main():
    import pandas as pd
    df = pd.DataFrame()
    for n in range(8, 98, 1):
        df = df.append(simulation(100, n))
        print(n)


if __name__ == '__main__':
    main()
    
