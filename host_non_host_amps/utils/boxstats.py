def return_clean(df):    
    ## cleaning off outliers
    import pandas as pd
    import matplotlib.pyplot as plt
    from .general_functions import remove_outlier

    mammal = remove_outlier(df[df.status == 'mammal'], 'normalized')
    mammal = mammal[['status', 'normalized']]
    nonmammal = remove_outlier(df[df.status == 'non-mammal'], 'normalized')
    nonmammal = nonmammal[['status', 'normalized']]
    plant = remove_outlier(df[df.status == 'plant'], 'normalized')
    plant = plant[['status', 'normalized']]
    environmental = remove_outlier(df[df.status == 'environmental'], 'normalized')
    environmental = environmental[['status', 'normalized']]

    pd.concat([mammal,
               nonmammal,
               plant,
               environmental]).boxplot(column='normalized',
                                       by='status',
                                       grid=False)
                                       
    plt.ylabel('AMPs per million of assembled base pairs')
    plt.tight_layout()
    plt.savefig('analysis/normamps_offoutlier_hostvsenv.svg')

    print(f'''Samples per set
    Mammal: {len(mammal)},
    Non-mammal: {len(nonmammal)},
    Plant: {len(plant)},
    Environmental: {len(environmental)}''')

    return [mammal, nonmammal, plant, environmental]


def statistics(mammal, nonmammal, plant, environmental):
    ## testing the statistics
    # uses two-sided tests for each possible pair of lists:
    import pandas as pd
    from scipy.stats import mannwhitneyu
        
    g = mannwhitneyu(mammal.normalized, nonmammal.normalized)
    print(f'Mammal vs. non mammal: {g}')
    g = mannwhitneyu(mammal.normalized, plant.normalized)
    print(f'Mammal vs. plant: {g}')
    g = mannwhitneyu(mammal.normalized, environmental.normalized)
    print(f'Mammal vs. environmental: {g}')
    g = mannwhitneyu(nonmammal.normalized, plant.normalized)
    print(f'Non mammal vs. plant: {g}')
    g = mannwhitneyu(nonmammal.normalized, environmental.normalized)
    print(f'Non mammal vs. environmental: {g}')
    g = mannwhitneyu(plant.normalized, environmental.normalized)
    print(f'Plant vs. environmental: {g}')


def plot_test(df):
    print('... plotting AMPs normalized richness')
    mammal, nonmammal, plant, environmental = return_clean(df)
    print('... testing statistic relevance')
    statistics(mammal, nonmammal, plant, environmental)

