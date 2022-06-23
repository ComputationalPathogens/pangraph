import xgboost as xgb
import pandas as pd
import numpy as np
import shap
import matplotlib.pyplot as plt
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import LabelEncoder
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.metrics import precision_recall_fscore_support, accuracy_score

def model_eval(predict, label):
    """
    Parameters
    ----------
    predict : Predictions made by model as a list
    label : Correct labels corresponding to predictions as a list
    Returns
    -------
    accuracy : Accuracy of model with supplied predictions

    """
    accuracy = np.sum([predict[i]==label[i] for i in range(len(predict))])/len(predict)

    return accuracy

def load_models(modelnums):
    models = []
    features = []
    labels = []
    for num in modelnums:
        bst = xgb.Booster()
        file = str(num) + '.model'
        filen = str(num) + '.npy'
        bst.load_model(file)
        temp = np.load(filen, allow_pickle=True)
        models.append(bst)
        features.append(temp.item().get('features'))
        labels.append(temp.item().get('labels'))

    return models, features, labels

def load_data(dataloc, filenamenp = '/processed_data/featuresfiltered.pkl', filenamecsv = '/processed_data/counts.csv'):
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
    datapth = dataloc + filenamenp
    labelpth = dataloc + filenamecsv
    #data = np.load(datapth, allow_pickle=True)
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


def train_model(k, features, labels, unencoded_labels, save, datadir):
    """
    k - amount of folds if doing cross fold validation (1 if not)
    features - x_train
    labels = y_train
    params - model parameters
    save - true to save models, false if not saving, also saves test data fold for accompanying model
    """
    params = {'objective':'multi:softmax', 'num_class': '11', 'max_depth': '12'}
    splits = np.load('/home/liam/compare/processed_data/foldsplits.npy', allow_pickle=True)
    #unknown = pd.read_pickle(datadir + '/processed_data/unknownfeatures.pkl')
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
        #utrain = sk_obj.transform(unknown)
        featmask = sk_obj.get_support()
        featnames = features.columns[featmask]
        xgb_matrix = xgb.DMatrix(Xtrain, label=Ytrain, feature_names=featnames)
        booster = xgb.train(params, xgb_matrix)
        xgb_test = xgb.DMatrix(Xtest, feature_names=featnames)
        #xgb_unk = xgb.DMatrix(utrain, feature_names=featnames)

        final_models.append(booster)
        final_features.append(xgb_test)
        final_labels.append(Ytest)
        final_train.append(Xtrain)
        final_train_y.append(Ytrain)
        #final_unknown.append(xgb_unk)
        shapimp = shap.TreeExplainer(booster)
        shap_vals = shapimp.shap_values(xgb_matrix)
        shap.summary_plot(shap_vals, Xtrain, show=False)
        plt.savefig('shapplot.png')
        #sortimp = {k: v for k, v in sorted(getimp.items(), key=lambda item: item[1], reverse=True)}

    test_model(final_models, final_features, final_labels, unencoded_labels, final_unknown, datadir)
    return final_models, final_features, final_labels

def test_model(final_models, final_features, final_labels, labels_unencoded, final_unknown, datadir):
    count = 0
    for model, xtest, ytest, utest in zip(final_models, final_features, final_labels, final_unknown):
        count+=1
        xgb_test_matrix = xtest
        prediction = model.predict(xgb_test_matrix)
        accuracy = accuracy_score(ytest, prediction)
        prec_recall = precision_recall_fscore_support(ytest, prediction, average=None)
        prec_recall = np.transpose(prec_recall)
        prec_recall = pd.DataFrame(data=prec_recall, index=labels_unencoded, columns=['Precision','Recall','F-Score','Supports'])
        model_report = datadir + '/processed_data/' + str(count) + 'summary.csv'
        prec_recall.to_csv(model_report)
        print(accuracy)
        print(prec_recall)
        print("UNKNOWN")
        #prediction = model.predict(utest)
        #with open(datadir + '/processed_data/fold' + str(count) + 'predicts.txt','w') as file:
        #    for x in range(len(prediction)):
        #        file.write(str(x) + '\t' + str(prediction[x]) + '\n')
        print(prediction)
    return
