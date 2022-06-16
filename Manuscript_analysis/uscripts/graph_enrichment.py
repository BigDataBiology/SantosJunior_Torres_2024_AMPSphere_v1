def main():
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    data = pd.read_table('data/adjust_significant_function.csv.xz', sep=',')
    data = data.sort_values(by='times')
    
    sns.barplot(data=data.tail(10), x='times', y='eggnog_OG', orient='h', color='black')
    plt.xlabel('Enrichment (fold)')
    plt.ylabel('EggNOG ortholog groups')
    plt.tight_layout()
    plt.savefig('top_og_enrichment.svg')
    plt.close()
    
    sns.barplot(data=data.tail(10), x='count_AMP', y='eggnog_OG', orient='h', color='black')
    plt.xlabel('AMP orthologs')
    plt.ylabel('EggNOG ortholog groups')
    plt.tight_layout()
    plt.savefig('most_frequent_og.svg')
    plt.close()
    
    data = pd.read_table('data/amp_COG.tsv.xz', names=['cog_class', 'count_AMP', 'fraction'])
    data = data.sort_values(by='count_AMP')
    
    sns.barplot(data=data.tail(10), x='count_AMP', y='cog_class', orient='h', color='black')
    plt.ylabel('COG class')
    plt.xlabel('AMP orthologs')
    plt.tight_layout()
    plt.savefig('most_frequent_cogclass.svg')


if __name__ == '__main__':
    main()
    
