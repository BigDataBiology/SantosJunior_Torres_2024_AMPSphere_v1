import pandas as pd
from itertools import chain
from utils.environment_classification import environments_list

from sklearn.model_selection import train_test_split
from sklearn.multiclass import OneVsRestClassifier

envdict = dict()
for ev, vv in zip(['non_mammalian', 'mammalian_gut', 'mammalian_other', 
                   'built_env.', 'freshwater', 'extreme',
                   'wastewater', 'soil', 'plant associated',
                   'marine', 'other'], environments_list):
        for i in vv: envdict[i] = ev


def gm(tp, fn, fp, tn):
    acc = (tp+tn)/(tp+fp+fn+tn)
    sens = tp/(tp+fn)
    sp = tn/(tn+fp)
    prec = tp/(tp+fp)
    f1 = 2*(sens*prec)/(sens+prec)
    n1 = (tp*tn) - (fp*fn)
    n2 = (tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)
    mcc = n1/n2
    return acc, sens, sp, prec, f1, mcc
    
    
def confusion_scores(y_test, predictions, labels):
    import pandas as pd
    from sklearn.metrics import confusion_matrix
    
    confusion = pd.DataFrame(confusion_matrix(y_test, predictions),
                             columns=labels.keys(),
                             index=labels.keys())
    
    ax1 = confusion.sum(axis=1)
    ax0 = confusion.sum(axis=0)
    
    classifier_performance = []
    for i in labels:
        tp = confusion.loc[i, i]
        fn = ax1.loc[i] - tp
        fp = ax0.loc[i] - tp
        tn = ax0.sum() - tp - fn - fp
        gmpred = gm(tp, fn, fp, tn)
        mat = [i]
        for j in [tp, fn, fp, tn]: mat.append(j)
        for j in gmpred: mat.append(j)
        classifier_performance.append(mat)
    
    df = pd.DataFrame(classifier_performance,
                      columns=['sample', 'tp', 'fn',
                               'fp', 'tn', 'accuracy',
                               'sensitivity', 'specificity',
                               'precision', 'f1-score',
                               'MCC'])
    return df


gen = pd.read_table('data/gmsc_amp_genes_envohr_source.tsv.gz')
gen = gen[['amp', 'general_envo_name']].drop_duplicates()
gen = gen[~(gen.general_envo_name.isna())]
gen.general_envo_name = [envdict[x] for x in gen.general_envo_name]
gen = gen.drop_duplicates()

motifs = pd.read_table('analysis/AMPSphere_v.2022-03.annotation.tsv')
motifs.rename({'id': 'amp'}, axis=1, inplace=True)
gen = gen.merge(on='amp', right=motifs)
gen = gen.drop(['amp', 'AA_rich', 'AA_absent', 'H', 'log2(gene_prob)'], axis=1)
gen = gen.dropna()
gen.motif_match = [x.split('|') for x in gen.motif_match]

motifs = set(chain.from_iterable(gen.motif_match))
for i in motifs:
    gen[i] = [1 if i in x else 0 for x in gen.motif_match]

gen = gen.drop('motif_match', axis=1)
gen = gen.replace([True, False], [1, 0])
gen.rename({'general_envo_name': 'label'}, axis=1, inplace=True)
gen = gen.sort_values(by='label')

Y = gen['label']
labels = dict()
for idx, val in enumerate(set(gen['label'])): labels[val] = idx
print(labels)

Y = [labels[x] for x in Y]
X = gen[gen.columns[1:]]
X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.3, random_state=42)

classifier = OneVsRestClassifier(RandomForestClassifier(n_estimators=501, criterion='entropy', random_state=42)).fit(X_train, y_train)
predictions = classifier.predict(X_test)
confusion_scores(y_test, predictions, labels)

