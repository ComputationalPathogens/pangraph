import pandas as pd
import numpy as np
import torch
import torch.nn.functional as F
from Bio import Seq, SeqIO
from torch.nn import Linear
from torch_geometric.nn import (GraphConv, SAGPooling, global_mean_pool, JumpingKnowledge)
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import precision_recall_fscore_support, accuracy_score
from captum.attr import Saliency, IntegratedGradients
from pyfiles import importance
from datetime import datetime

device = torch.device('cuda')

####
# Model definition for graph network
####

class SAGPool(torch.nn.Module):
    def __init__(self, num_layers, hidden, num_class, ratio=0.1):
        super(SAGPool, self).__init__()
        self.conv1 = GraphConv(656, hidden, aggr='mean')
        self.convs = torch.nn.ModuleList()
        self.pools = torch.nn.ModuleList()
        self.convs.extend([
            GraphConv(hidden, hidden, aggr='mean')
            for i in range(num_layers - 1)
        ])
        self.pools.extend(
            [SAGPooling(hidden, ratio) for i in range((num_layers) // 2)])
        self.jump = JumpingKnowledge(mode='cat')
        self.lin1 = Linear(num_layers * hidden, hidden)
        self.lin2 = Linear(hidden, num_class)

    def reset_parameters(self):
        self.conv1.reset_parameters()
        for conv in self.convs:
            conv.reset_parameters()
        for pool in self.pools:
            pool.reset_parameters()
        self.lin1.reset_parameters()
        self.lin2.reset_parameters()

    def forward(self, data):
        x, edge_index, batch = data.x, data.edge_index, data.batch
        x = F.relu(self.conv1(x, edge_index))
        xs = [global_mean_pool(x, batch)]
        for i, conv in enumerate(self.convs):
            x = F.relu(conv(x, edge_index))
            xs += [global_mean_pool(x, batch)]
            if i % 2 == 0 and i < len(self.convs) - 1:
                pool = self.pools[i // 2]
                x, edge_index, _, batch, _, _ = pool(x, edge_index,
                                                     batch=batch)
        x = self.jump(xs)
        x = F.relu(self.lin1(x))
        x = F.dropout(x, p=0.5, training=self.training)
        x = self.lin2(x)
        return F.log_softmax(x, dim=-1)

    def __repr__(self):
        return self.__class__.__name__

def get_acc(loader, model):
    model.eval()
    ytrue, ypred = [],[]
    for d in loader:
        ytrue.append(d.y[0].item())
        d = d.to(device)
        m = model.to(device)
        pred = m(d)
        pred = pred.max(dim=1)[1]
        ypred.append(pred[0].item())
    accuracy = accuracy_score(ytrue, ypred)
    return accuracy
    

def build(datadir, meta, metapth):
    ####
    # Turn into helper function??
    ####
    meta = False
    splits = np.load(datadir + '/processed_data/foldsplits.npy', allow_pickle=True)
    colnames = ['id', 'assembly', 'genus', 'species', 'seqfile', 'cntfile', 'meta']
    samples = pd.read_csv(datadir + '/processed_data/clean.csv', names=colnames)
    species = samples.species.tolist()
    enc = 0
    indspec = {}
    specdict = {}
    labels_unencoded = []
    #TURN INTO HELPER FUNCTION IN DIFFERENT CLASS EVERYWHERE THIS EXISTS
    for s in set(species):
        labels_unencoded.append(s)
    labels_unencoded.sort()
    for s in labels_unencoded:
        specdict[s] = enc
        indspec[enc] = s
        enc += 1
    
    #labels_unencoded.append('rev_mel')
    predictdict = {}
    graphpredicts = {}
    print("NUMCLASSES:%i" % (enc))
    print("Specdict")
    print(specdict)
    print(labels_unencoded)
    ####
    # Just printing general information like which number represents which class
    ####
    final_models = []
    final_features = []
    final_labels = []
    final_train = []
    final_train_y = []
    final_unknown = []
    final_meta = []
    final_importance = [[],[],[],[],[]]
    
    for fold in range(5):
        print("Fold#: ", str(fold), datetime.now())
        donegraphs = torch.load(datadir + '/processed_data/fold' + str(fold+1) + 'dataset.pkl')
        if meta:
            metagraphs = torch.load(datadir + '/processed_data/' + metapth + 'fold' + str(fold+1) + 'dataset.pkl')
        #unknowngraphs = torch.load(datadir + '/processed_data/unknown/fold' + str(fold+1) + 'dataset.pkl')
        model = SAGPool(3, 512, enc)
        train_data = []
        test_data = []
        meli = []
        abor = []
        suis = []
        cani = []
        ovis = []
        importance_graphs = []
        
        ####
        # Splitting fold dataset into train/test sets based on saved splits
        ####
        for g in donegraphs:
            if int(g.graphind) in splits[0][fold]:
                train_data.append(g)
            else:
                test_data.append(g)
                
        """
        for x in range(len(metagraphs)):
            if x % 2 == 0:
                train_data.append(metagraphs[x])
            else:
                test_data.append(metagraphs[x])
        print(x)
        """    
        ####
        # Seperating graphs we want feature importance for
        ####
        for g in donegraphs:
            if g.y.item() == 5:
                meli.append(g)
            elif g.y.item() == 0:
                abor.append(g)
            elif g.y.item() == 8:
                suis.append(g)
            elif g.y.item() == 2:
                cani.append(g)
            elif g.y.item() == 6:
                ovis.append(g)
        importance_graphs = [meli, abor, suis, cani, ovis]
        
        ####
        # Training each folds model for 50 epochs, outputting accuracy of train/test set at each iteration
        ####
        train_loader = DataLoader(train_data, batch_size=1, shuffle=True)
        test_loader = DataLoader(test_data, batch_size=1, shuffle=True)
        #unknown_loader = DataLoader(unknowngraphs, batch_size=1)
        if meta:
            meta_loader = DataLoader(metagraphs, batch_size=1, shuffle=False)
        device = torch.device('cuda')
        model = model.to(device)
        optimizer = torch.optim.Adam(model.parameters(), lr=0.0001)
        for epoch in range(150):
            model.train()
            loss_all = 0
            for d in train_loader:
                d = d.to(device)
                optimizer.zero_grad()
                output = model(d)
                loss = F.nll_loss(output, d.y)
                loss.backward()
                loss_all += loss.item()
                optimizer.step()
            train_acc = get_acc(train_loader, model)
            test_acc = get_acc(test_loader, model)
            print("EPOCH #%i  FOLD#%i  LOSS:%i  TRAINACC:%f  TESTACC:%f" % (epoch+1,fold+1, loss_all, train_acc, test_acc))
        
        ####
        # After this folds model is trained, pass it to feature importance module to extract important class features
        ####

        specnames = ['melitensis', 'abortus', 'suis', 'canis', 'ovis']
        #for ig in range(len(importance_graphs)):
            #print(specnames[ig])
        #final_importance[1].append(importance.importance(datadir, importance_graphs[1], model, int(fold),specnames[1]))
 
        final_models.append(model)
        final_features.append(test_loader)
        final_train.append(train_loader)
        #final_unknown.append(unknown_loader)
        if meta:
            final_meta.append(meta_loader)
        else:
            final_meta.append('')
        
         
    ####
    # Once training is finished we write important features to file and run testing splits on each model
    ####
    """
    print("TRAINDONE:", datetime.now())
    for x in range(5):
        with open(datadir + '/processed_data/' + str(x) + '_features.txt', 'w') as f:
            for feat in final_importance[x]:
                f.write(feat + '\n')
    """
    print(specdict)
    print("TRAINING DONE")
    ind = 0
    allmpred = []
    for m, xtest, mtest in zip(final_models, final_features, final_meta):
        print("TESTING: ", datetime.now())
        wrongdict = {}
        print("TESTING FOLD#%i" % (ind+1))
        ind +=1
        model.eval()
        correct = 0
        ytrue = []
        ypred = []
        #upred = []
        mpred = []
        if meta:
            for d in mtest:
                print("SAMPLE")
                d = d.to(device)
                m = m.to(device)
                pred = m(d)
                print(pred)
                mpred.append(pred)
                predmax = pred.max(dim=1)[1]
                print("METAPREDICTION: ", d.graphind, predmax)
            allmpred.append(mpred)

        ####
        # Apply models to unknown dataset and record predictions in file for comparision
        ####
        """
        with open(datadir + '/processed_data/fold' + str(ind) + 'predicts.txt', 'w') as file:
            for d in utest:
                d = d.to(device)
                m = m.to(device)
                pred = m(d)
                pred = pred.max(dim=1)[1]
                upred.append(pred[0].item())
                file.write(str(d.graphind.item()) + '\t' + str(pred[0].item()) + '\n')
        print("PREDICTIONS ON UKNNOWN SET")
        print(upred)
        """
        ####
        # Applying model to regular testing set
        ####
        for d in xtest:
            ytrue.append(d.y[0].item())
            d = d.to(device)
            m = m.to(device)
            pred = m(d)
            pred = pred.max(dim=1)[1]
            ypred.append(pred[0].item())
            if pred[0].item() not in predictdict:
                predictdict[pred[0].item()] = 1
            else:
                predictdict[pred[0].item()] += 1
            if d.graphind not in graphpredicts:
                graphpredicts[d.graphind.item()] = [pred[0].item()]
            else:
                graphpredicts[d.graphind.item()].append(pred[0].item())
            if pred[0].item() != d.y[0].item():
                if d.y[0].item() not in wrongdict:
                    wrongdict[d.y[0].item()] = 1
                else:
                    wrongdict[d.y[0].item()] += 1
            correct += pred.eq(d.y).sum().item()
        accuracy = accuracy_score(ytrue, ypred)
        prec_recall = precision_recall_fscore_support(ytrue,ypred)
        prec_recall = np.transpose(prec_recall)
        prec_recall = pd.DataFrame(data=prec_recall, index=labels_unencoded, columns=['Precision','Recall','F-Score','Supports'])
        model_report = datadir + '/processed_data/' + str(ind) + 'summary.csv'
        print(model_report)
        print(prec_recall)
        prec_recall.to_csv(model_report)
        print("Accuracy: %f" % (accuracy))
        print("Correct:#%i" %  correct)
        print("Total#%i" % len(xtest.dataset))
        
    if meta:
        print("METAGENOME BIN PREDICTIONS:")
        tpred = []
        for x in range(len(allmpred[0])):
            tpred.append(allmpred[0][x])
        tmax = []
        for x in range(4):
            for y in range(len(tpred)):
                tpred[y] += allmpred[x+1][y]
        print(tpred)
        for x in range(len(tpred)):
            tmax.append(tpred[x].max(dim=1)[1])
        print(tmax)
    print(predictdict)
