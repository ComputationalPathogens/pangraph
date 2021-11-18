import pandas as pd
import numpy as np
import torch
import torch.nn.functional as F
from torch.nn import Linear
from torch_geometric.nn import (GraphConv, SAGPooling, global_mean_pool, JumpingKnowledge)
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import precision_recall_fscore_support, accuracy_score
from captum.attr import Saliency, IntegratedGradients

device = torch.device('cuda')

class SAGPool(torch.nn.Module):
    def __init__(self, num_layers, hidden, num_class, ratio=0.8):
        super(SAGPool, self).__init__()
        self.conv1 = GraphConv(798, hidden, aggr='mean')
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

def model_forward(node_mask, data, model):
    integrated = Data(x=node_mask, edge_index=data.edge_index, y=data.y, batch=data.batch, ptr = data.ptr, graphind = data.graphind)
    out = model(integrated)
    return out

def explain(data, model, target=0):
    input_mask = torch.ones((data.x.shape[0],data.x.shape[1])).to(device)
    input_mask = torch.tensor(input_mask)
    ig = IntegratedGradients(model_forward)
    mask = ig.attribute(input_mask, target=data.y[0],additional_forward_args=(data,model,),internal_batch_size=1)
    node_mask = np.abs(mask.cpu().detach().numpy())
    if node_mask.max() > 0:
        node_mask = node_mask / node_mask.max()
    return node_mask
    
def build(datadir):
    splits = np.load(datadir + '/processed_data/foldsplits.npy', allow_pickle=True)
    colnames = ['id', 'assembly', 'genus', 'species', 'seqfile', 'cntfile', 'meta']
    samples = pd.read_csv(datadir + '/processed_data/clean.csv', names=colnames)
    species = samples.species.tolist()
    enc = 0
    indspec = {}
    specdict = {}
    labels_unencoded = []
    for s in set(species):
        labels_unencoded.append(s)
    labels_unencoded.sort()
    for s in labels_unencoded:
        specdict[s] = enc
        indspec[enc] = s
        enc += 1
    print("NUMCLASSES:%i" % (enc))
    print("Specdict")
    print(specdict)
    print(labels_unencoded)
    final_models = []
    final_features = []
    final_labels = []
    final_train = []
    final_train_y = []
    for fold in range(5):
        class_dist = {}
        donegraphs = torch.load(datadir + '/processed_data/fold' + str(fold+1) + 'dataset.pkl')
        model = SAGPool(6, 512, enc)
        train_data = []
        test_data = []
        meli = []
        presmatrix = []
        meligraphinds = []
        for ind in range(len(splits[1][fold])):
            if specdict[species[splits[1][fold][ind]]] == 0:
                meligraphinds.append(ind)
        for g in donegraphs:
            if g.y.item() not in class_dist:
                class_dist[g.y.item()] = 1
            else:
                class_dist[g.y.item()] += 1
            if int(g.graphind) in splits[0][fold]:
                train_data.append(g)
            else:
                if g.y.item() == 0:
                    meli.append(g)
                test_data.append(g)
        ind = 0
        for g in meli:
            pos,neg = 0,0
            presmatrix.append([])
            for node in range(len(g.x)):
                for posind in meligraphinds:
                    if g.x[node][posind].item() == 0:
                        neg += 1
                    else:
                        pos += 1
                if pos >= neg:
                    presmatrix[ind].append(1)
                else:
                    presmatrix[ind].append(0)
            g.meli = presmatrix[ind]
            #print(presmatrix[ind])
            ind += 1
        print(class_dist)
        train_loader = DataLoader(train_data, batch_size=1)
        test_loader = DataLoader(test_data, batch_size=1)
        expl_loader = DataLoader(meli, batch_size=1)
        device = torch.device('cuda')
        model = model.to(device)
        optimizer = torch.optim.Adam(model.parameters(), lr=0.0001)
        for epoch in range(100):
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
            print("EPOCH #%i  FOLD#%i  LOSS:%i" % (epoch+1,fold+1, loss_all))
            
            
        model.eval()
        for d in expl_loader:
            d = d.to(device)
            expl = explain(d, model)
            importance = np.zeros(expl.shape[0])
            maximp = 0
            maxind = 0
            for node in range(len(expl)):
                importance[node] = np.sum(expl[node])
                if importance[node] > maximp:
                    maximp = importance[node]
                    maxind = node
            #print(d.graphind)
            #print(d.y)
            count = 0
            seqs = []
            with open(datadir + '/processed_data/fasta/graph' + str(d.graphind.item()) + '.fasta') as f:
                for l in f:
                    if count % 2 != 0:
                        seqs.append(str.rsplit(l)[0])
                    count += 1
         
            print("HIGHEST IMPORTANCE SEQUENCE")
            print(maxind)
            print(d.meli[maxind].item())
            print(seqs[maxind])
                
        final_models.append(model)

        final_features.append(test_loader)
        final_train.append(train_loader)

    print(specdict)
    print("TRAINING DONE")
    ind = 0
    for m, xtest in zip(final_models, final_features):
        
        wrongdict = {}
        print("TESTING FOLD#%i" % (ind+1))
        ind +=1
        model.eval()
        correct = 0
        ytrue = []
        ypred = []
        for d in xtest:
            ytrue.append(d.y[0].item())
            d = d.to(device)
            m = m.to(device)
            pred = m(d)
            pred = pred.max(dim=1)[1]
            ypred.append(pred[0].item())
            if pred[0].item() != d.y[0].item():
                if d.y[0].item() not in wrongdict:
                    wrongdict[d.y[0].item()] = 1
                else:
                    wrongdict[d.y[0].item()] += 1
            #print("HIGHEST:%i vs. ACTUAL:%i" % (pred,d.y[0]))
            correct += pred.eq(d.y).sum().item()
        accuracy = accuracy_score(ytrue, ypred)
        prec_recall = precision_recall_fscore_support(ytrue,ypred)
        prec_recall = np.transpose(prec_recall)
        prec_recall = pd.DataFrame(data=prec_recall, index=labels_unencoded, columns=['Precision','Recall','F-Score','Supports'])
        model_report = datadir + '/processed_data/' + str(ind) + 'summary.csv'
        print(model_report)
        prec_recall.to_csv(model_report)
        print(prec_recall)
        print("Accuracy: %f" % (accuracy))
        print("Correct:#%i" %  correct)
        print("Total#%i" % len(xtest.dataset))