def species_amps_distribution(data):    
    ## creating species distribution over environments
    import pandas as pd
    from .general_functions import mergefixed

    data = data[['amp', 'fixed', 'general_envo_name']]
    data = data.drop_duplicates()

    new = dict()
    for record in data.groupby('general_envo_name'):
        new[record[0]] = mergefixed(record[1].fixed)

    species_distribution = pd.DataFrame.from_dict(new,
                                                  orient='columns')
    species_distribution = species_distribution.fillna(0)
    species_distribution.to_csv('analysis/species_distribution_by_env.tsv',
                                sep='\t',
                                header=True,
                                index=True)


def diversity_per_sample(data, spheres):
    ## calculating genera diversity per sample
    import pandas as pd
    from collections import Counter
    from .general_functions import mergefixed
    from skbio.diversity.alpha import shannon
    
    data = data.merge(on='amp', right=spheres)
    new = []
    for i in data.groupby(['general_envo_name', 'sample']):
        diversity_sp = mergefixed(i[1].fixed)
        diversity_sp = pd.DataFrame.from_dict(diversity_sp,
                                              orient='index')
        diversity_sp = diversity_sp[0]
        diversity_sp = shannon(diversity_sp)
        diversity_amp = Counter(i[1].family)
        diversity_amp = pd.DataFrame.from_dict(diversity_amp,
                                               orient='index')
        diversity_amp = diversity_amp[0]
        diversity_amp = shannon(diversity_amp)
        new.append([i[0][0],
                    i[0][1],
                    diversity_amp, diversity_sp])

    divshannon = pd.DataFrame(new,
                              columns=['status',
                                       'sample',
                                       'diversity_amp',
                                       'diversity_sp'])
                                       
    divshannon.to_csv('analysis/shannon_diversity_per_sample.tsv', sep='\t', header=True, index=None)


def distributions(data, spheres):
    print('... checking species distribution')
    species_amps_distribution(data)
    print('... checking diversity of amps')
    diversity_per_sample(data, spheres)

