import xgboost as xgb
import pandas as pd
import numpy as np
import shap
import matplotlib.pyplot as plt
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import LabelEncoder
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.metrics import precision_recall_fscore_support, accuracy_score

def load_data(dataloc, dname):
    """
    Parameters
    ----------
    dataloc : Base directory of nextflow execution

    Returns
    -------
    data : Feature matrix of [NumSamples]*[NumFeatures] shape
    labels_encoded : Labels corresponding to feature matrix [Numsamples] length
    label_encoder.classes_ : the unencoded classes the model is being trained on
    
    """
    datapth = dataloc + '/processed_data/' + str(dname) + '_featuresfiltered.pkl'
    labelpth = dataloc + '/processed_data/' + str(dname) + '_counts.csv'
    data = pd.read_pickle(datapth)
    colnames = ['id', 'assembly', 'genus', 'species', 'seqfile', 'cntfile', 'meta']
    labels = pd.read_csv(labelpth, names=colnames)
    labels = labels.species.tolist()
    specdict = {}
    encoded = []
    labels_unencoded = []
    enc = 0
    for s in set(labels):
        labels_unencoded.append(s)
    labels_unencoded.sort()
    for s in labels_unencoded:
        specdict[s] = enc
        enc += 1
    for l in labels:
        encoded.append(specdict[l])
        
    return data, encoded, labels_unencoded


def train_model(k, features, labels, unencoded_labels, save, datadir, dname):
    """
    k - amount of folds if doing cross fold validation (1 if not)
    features - x_train
    labels = y_train
    params - model parameters
    save - true to save models, false if not saving, also saves test data fold for accompanying model
    """
    params = {'objective':'multi:softmax', 'num_class': '11', 'max_depth': '12'}
    splits = np.load(datadir + '/processed_data/' + str(dname) + '_foldsplits.npy', allow_pickle=True)
    count = 0
    num_feats = 2000000
    final_models = []
    final_features = []
    final_labels = []
    final_train = []
    final_train_y = []
    final_unknown = []
    for x in range(5):
        count+=1
        combinesplits = [*splits[0][x], *splits[1][x]]
        testsplits = [*splits[2][x]]
        sk_obj = SelectKBest(f_classif, k=num_feats)
        print(combinesplits)
        print(testsplits)
        print(features)
        Xtrain,Xtest = features.iloc[combinesplits], features.iloc[testsplits]
        Ytrain = [labels[i] for i in combinesplits]
        Ytest = [labels[i] for i in testsplits]

        Xtrain = sk_obj.fit_transform(Xtrain, Ytrain)
        Xtest = sk_obj.transform(Xtest)
        featmask = sk_obj.get_support()
        featnames = features.columns[featmask]
        xgb_matrix = xgb.DMatrix(Xtrain, label=Ytrain, feature_names=featnames)
        booster = xgb.train(params, xgb_matrix)
        xgb_test = xgb.DMatrix(Xtest, feature_names=featnames)

        final_models.append(booster)
        final_features.append(xgb_test)
        final_labels.append(Ytest)
        final_train.append(Xtrain)
        final_train_y.append(Ytrain)
        """
        #Feature Importance Code
        shapimp = shap.TreeExplainer(booster)
        shap_vals = shapimp.shap_values(xgb_matrix)
        shap.summary_plot(shap_vals, Xtrain, show=False)
        plt.savefig('shapplot.png')
        sortimp = {k: v for k, v in sorted(getimp.items(), key=lambda item: item[1], reverse=True)}
        """

    test_model(final_models, final_features, final_labels, unencoded_labels, datadir, dname)
    return final_models, final_features, final_labels

def test_model(final_models, final_features, final_labels, labels_unencoded, datadir, dname):
    count = 0
    for model, xtest, ytest in zip(final_models, final_features, final_labels):
        count+=1
        xgb_test_matrix = xtest
        prediction = model.predict(xgb_test_matrix)
        accuracy = accuracy_score(ytest, prediction)
        prec_recall = precision_recall_fscore_support(ytest, prediction, average=None)
        prec_recall = np.transpose(prec_recall)
        prec_recall = pd.DataFrame(data=prec_recall, index=labels_unencoded, columns=['Precision','Recall','F-Score','Supports'])
        model_report = datadir + '/processed_data/' + str(dname) + '_' + str(count) + 'summary.csv'
        prec_recall.to_csv(model_report)
        print(accuracy)
        print(prec_recall)
    return
