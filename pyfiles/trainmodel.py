import xgboost as xgb
import pandas as pd
import numpy as np
from tensorflow import keras
from tensorflow.keras.layers import Dense, Dropout
from keras_tuner import Hyperband
#import tensorflow as tf
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

def load_data(dataloc, filenamenp = '/processed_data/features.pkl', filenamecsv = '/processed_data/cleanwcounts.csv'):
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
    count = 0
    num_feats = 2000000
    final_models = []
    final_features = []
    final_labels = []
    final_train = []
    final_train_y = []
    for x in range(5):
        print(features.shape)
        count+=1
        sk_obj = SelectKBest(f_classif, k=num_feats)
        Xtrain,Xtest = features.iloc[splits[0][x]], features.iloc[splits[2][x]]
        Ytrain = [labels[i] for i in splits[0][x]]
        Ytest = [labels[i] for i in splits[2][x]]

        Xtrain = sk_obj.fit_transform(Xtrain, Ytrain)
        Xtest = sk_obj.transform(Xtest)
        featmask = sk_obj.get_support()
        featnames = features.columns[featmask]
        xgb_matrix = xgb.DMatrix(Xtrain, label=Ytrain)
        booster = xgb.train(params, xgb_matrix)
        xgb_test = xgb.DMatrix(Xtest)

        final_models.append(booster)
        final_features.append(Xtest)
        final_labels.append(Ytest)
        final_train.append(Xtrain)
        final_train_y.append(Ytrain)
    test_model(final_models, final_features, final_labels, unencoded_labels, datadir)
    return final_models, final_features, final_labels

def test_model(final_models, final_features, final_labels, labels_unencoded, datadir):
    count = 0
    for model, xtest, ytest in zip(final_models, final_features, final_labels):
        count+=1
        xgb_test_matrix = xgb.DMatrix(xtest)
        prediction = model.predict(xgb_test_matrix)
        accuracy = accuracy_score(ytest, prediction)
        prec_recall = precision_recall_fscore_support(ytest, prediction, average=None)
        prec_recall = np.transpose(prec_recall)
        prec_recall = pd.DataFrame(data=prec_recall, index=labels_unencoded, columns=['Precision','Recall','F-Score','Supports'])
        model_report = datadir + '/processed_data/' + str(count) + 'summary.csv'
        prec_recall.to_csv(model_report)
        print(accuracy)
        print(prec_recall)
    return

def build_model(hp):
    model = keras.Sequential()
    model.add(
        Dense(
             units=hp.Int("units", min_value=32, max_value=512, step=32),
            activation="relu", input_dim=1000
        )
    )
    model.add(Dropout(0.5))
    model.add(Dense(14, kernel_initializer='uniform', activation='softmax'))
    model.compile(
        optimizer=keras.optimizers.Adam(
            hp.Choice("learning_rate", values=[1e-2, 1e-3, 1e-4])
        ),
        loss="sparse_categorical_crossentropy",
        metrics=["accuracy"],
    )

    return model


def train_keras(k, data, label_encoded_y, labels_unencoded):
    num_feats = 1000
    kf = StratifiedKFold(n_splits=5, shuffle=True)
    final_models = []
    final_features = []
    final_labels = []
    final_hps = []

    stop_early = keras.callbacks.EarlyStopping(monitor='val_loss', patience=6)
    fold = 0
    for train_index, test_index in kf.split(data, label_encoded_y):
        tuner = Hyperband(
        build_model,
        objective="val_accuracy",
        max_epochs=10,
        factor=3,
        overwrite=True,
        directory="/home/liam/pathogen-pipeline/hyp",
        project_name="BCP",
    )
        fold += 1
        Xtrain = data[train_index]
        Xtest = data[test_index]
        Ytrain = label_encoded_y[train_index]
        Ytest = label_encoded_y[test_index]
        sk_obj = SelectKBest(f_classif, k=num_feats)
        Xtrain = sk_obj.fit_transform(Xtrain, Ytrain)
        Xtest = sk_obj.transform(Xtest)
        tuner.search(Xtrain, Ytrain, epochs=20, validation_split=0.25,callbacks=[stop_early], verbose=0)
        best_hps = tuner.get_best_hyperparameters(num_trials=1)[0]
        model = tuner.hypermodel.build(best_hps)
        final_hps.append(best_hps)
        final_features.append(Xtest)
        final_labels.append(Ytest)

        history = model.fit(Xtrain, Ytrain, epochs=50, validation_split=0.25, verbose = 0)
        val_acc_per_epoch = history.history['val_accuracy']
        best_epoch = val_acc_per_epoch.index(max(val_acc_per_epoch)) + 1
        hypermodel = tuner.hypermodel.build(best_hps)
        hypermodel.fit(Xtrain, Ytrain, epochs=best_epoch, validation_split=0.25, verbose = 0)
        final_models.append(hypermodel)

    return final_hps, final_models, final_features, final_labels

def test_keras(final_models, final_features, final_labels):
    highest_acc = 0
    best_model = []
    for model, Xtest, Ytest in zip(final_models, final_features, final_labels):
        eval_result = model.evaluate(Xtest, Ytest, verbose = 0)
        print("[test loss, test accuracy]:", eval_result)
        if eval_result[1] > highest_acc:
            highest_acc = eval_result[1]
            best_model = model

    return best_model
