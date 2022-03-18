## stating functions
def aps(df, name):
    # AMPs per sample function
    # Takes as input:
    # df - a dataframe containing AMPs, sample, and general environment name
    # envolist - a list of environments in the general environment name
    # name - a string to name the dataset
    # returns a dataframe with the number of AMPs per sample
    # and a tag for the environment
    import pandas as pd
    from collections import Counter
    
    envolist = df[df['general_envo_name'] == name]
    envolist = envolist['sample']
    envolist = Counter(envolist)
    envolist = pd.DataFrame.from_dict(envolist,
                                      orient='index')
    envolist = envolist.reset_index()
    envolist = envolist.rename({'index': 'sample',
                                0: 'AMPs'},
                               axis=1)
    envolist['status'] = name
    return envolist


def remove_outlier(df_in, col_name):
    # Remove outliers function
    # Takes as input:
    # df - a dataframe containing normalized AMPs per sample and a tag ('status')
    # col_name - string naming the column of interest containing the values to be filtered
    # returns a dataframe with the rows with a value in the column of
    # interest inside the low and high fences
    import pandas as pd
    
    q1 = df_in[col_name].quantile(0.25)
    q3 = df_in[col_name].quantile(0.75)
    iqr = q3-q1 #Interquartile range
    fence_low  = q1-1.5*iqr
    fence_high = q3+1.5*iqr
    df_out = df_in.loc[(df_in[col_name] > fence_low) & (df_in[col_name] < fence_high)]
    return df_out


def mergefixed(list_of_str):
    # Merge sets given as strings
    # input is a list of set elements converted into string
    # return a dictionary-like object with the elements of the original stringed-set
    from collections import Counter

    final = []
    for i in list_of_str:
        i = i.replace('{', '').replace('}', '')
        i = i.replace("'", '').split(', ')
        for ij in i:
            final += i

    return Counter(final)

