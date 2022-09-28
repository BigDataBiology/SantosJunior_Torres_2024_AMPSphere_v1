def maplabels():
    import pandas as pd
    import seaborn as sns
    from matplotlib.colors import ListedColormap

    d1 = pd.read_table('merged_to_donor_samples.tsv')
    d2 = pd.read_table('merged_to_recipient_samples_post.tsv')
    d3 = pd.read_table('merged_to_recipient_samples_pre.tsv')

    s = set(d1.outcome)
    s = s.union(set(d2.outcome))
    s = s.union(set(d3.outcome))

    paired = sns.color_palette("Paired")
    cmap = dict()
    for n, i in enumerate(s): cmap[i] = paired[n]
    return cmap

def plot_per_genome(fin, normalize=None, name=None):
    import pandas as pd
    from matplotlib import pyplot as plt
    from collections import Counter
    if normalize == None: normalize = True
    if name == None: name = fin.split('/')[-1].replace('.tsv', '.svg')
    cmap = maplabels()
    data = pd.read_table(fin)
    df = data.groupby('#_amps')['outcome'].apply(lambda x: Counter(x))
    df = df.reset_index().pivot(columns='level_1', index='#_amps', values='outcome')
    df = df.fillna(0)
    if normalize:
        df = df.T * 100 / df.sum(axis=1)
        df.T.plot.bar(stacked=True, color=[cmap[c] for c in df.T.columns])
        plt.ylabel('Proportion of outcomes')
    else:
        df.plot.bar(stacked=True, color=[cmap[c] for c in df.columns])
        plt.ylabel('Metagenome assembled genomes')
    plt.title(name)
    plt.xlabel('AMPs per genome')
    plt.tight_layout()
    plt.savefig(name)
    plt.close()


def plot_genomes():
    plot_per_genome('merged_to_donor_samples.tsv', True, 'merged_to_donor_samples_norm.svg')
    plot_per_genome('merged_to_recipient_samples_post.tsv', True, 'merged_to_recipient_samples_post_norm.svg')
    plot_per_genome('merged_to_recipient_samples_pre.tsv', True, 'merged_to_recipient_samples_pre_norm.svg')
    plot_per_genome('merged_to_donor_samples.tsv', False, 'merged_to_donor_samples.svg')
    plot_per_genome('merged_to_recipient_samples_post.tsv', False, 'merged_to_recipient_samples_post.svg')
    plot_per_genome('merged_to_recipient_samples_pre.tsv', False, 'merged_to_recipient_samples_pre.svg')

