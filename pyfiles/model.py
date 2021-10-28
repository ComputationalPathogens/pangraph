import pandas as pd
import numpy as np
import random
import torch
import torch.nn.functional as F
from torch.nn import Linear
from torch_geometric.nn import (GraphConv, SAGPooling, global_mean_pool, JumpingKnowledge)
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader
from sklearn.model_selection import StratifiedKFold

save = True

class SAGPool(torch.nn.Module):
    def __init__(self, num_layers, hidden, ratio=0.8):
        super(SAGPool, self).__init__()
        self.conv1 = GraphConv(1001, hidden, aggr='mean')
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
        self.lin2 = Linear(hidden, 31)

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
    
def build(datadir):
    donegraphs = torch.load('/home/liam/pytorchtest/graphdataset.pkl')
    labels = []
    labelcounts = {}
    for graph in donegraphs:
        if graph.y.item() not in labelcounts:
            labelcounts[graph.y.item()] = 1
        else:
            labelcounts[graph.y.item()] += 1
    cleangraphs = []
    for graph in donegraphs:
        if labelcounts[graph.y.item()] >= 5:
            cleangraphs.append(graph)
            labels.append(graph.y.item())
    kf = StratifiedKFold(n_splits=5, shuffle=True)
    print(labelcounts)
    final_models = []
    final_features = []
    final_labels = []
    final_train = []
    final_train_y = []
    count = 0
    for train_index, test_index in kf.split(cleangraphs, labels):
        count += 1
        model = SAGPool(5, 256)
        train_data = [cleangraphs[i] for i in train_index]
        #train_y = labels[train_index]
        test_data = [cleangraphs[i] for i in test_index]
        #test_y = labels[test_index]
        train_loader = DataLoader(train_data, batch_size=1)
        test_loader = DataLoader(test_data, batch_size=1)
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
            print(loss_all)
        final_models.append(model)
        #if save == True:
        #    modelsave = '/home/liam/pytorchtest/GraphModel/' + str(count) + '.pkl'
        #    torch.save(modelsave, model)
        final_features.append(test_loader)
        #final_labels.append(test_y)
        final_train.append(train_loader)
        #final_train_y.append(train_y)
        
    for m, xtest in zip(final_models, final_features):
        model.eval()
        correct = 0
        print(xtest)
        for d in xtest:
            d = d.to(device)
            m = m.to(device)
            pred = m(d).max(dim=1)[1]
            #print(indspec[pred.item()], " vs. ", indspec[d.y[0].item()])
            correct += pred.eq(d.y).sum().item()
    
        print("Accuracy: %f" % (correct/len(xtest.dataset)))
        print((correct/len(xtest.dataset)))
        print("Correct:#%i" %  correct)
        print("Total#%i" % len(xtest.dataset))