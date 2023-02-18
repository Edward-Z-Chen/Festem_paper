#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 27 09:45:58 2022
This script is modified from the tutorial of the TN_test paper.
https://github.com/jessemzhang/tn_test/blob/master/tntest_tutorial.ipynb

@author: Edward Chen
"""
import time
import pickle
import itertools
# import pyreadr
import numpy as np
import pandas as pd
import scanpy as sc
from truncated_normal import truncated_normal as tn
from scipy.stats import ttest_ind
from pyreadr import read_r

# Load data
adata = read_r("./results/ifnb_counts_forpython.RData")
counts = adata["counts"]
features = counts.index.values
counts = counts.transpose()
counts = counts.to_numpy()
counts = counts.astype(int)
del adata

# Louvain
def scanpy_cluster(X, features, tsne=None, plot=False, resolution=1.0):
    "Scanpy implementation of Seurat's pipeline"

    adata = sc.AnnData(X=X, dtype="int")
    adata.var['genes_ids'] = features

    
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3")
    adata = adata[:, adata.var['highly_variable']]
    sc.pp.scale(adata)

    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=20, n_pcs=10)
    sc.tl.louvain(adata, resolution=resolution)

    if plot:
        if tsne is not None:
            adata.obsm['X_tsne'] = tsne
        else:
            sc.tl.tsne(adata)
        sc.pl.tsne(adata, color=['louvain'])

    labels = np.array(adata.obs['louvain'].astype(int))

    return labels


# Split the dataset
np.random.seed(0)
n = counts.shape[0]
inds1 = np.sort(np.random.choice(range(n), n // 2, replace=False))
inds2 = np.ones(n).astype(bool)
inds2[inds1] = False
X1, X2 = counts[inds1], counts[inds2]

samp_labels = np.array(['Partition 1' if i else 'Partition 2' for i in inds2])
# plot_labels_legend(tsne[:, 0], tsne[:, 1], samp_labels)

# Generate clusters using X1
labels1 = scanpy_cluster(X1, features, plot=False, resolution=1.11)

# Fit hyperplanes using X1 (one-versus-one)
# from sklearn.svm import SVC
# svm = SVC(kernel='linear', C=100)

# Fit hyperplanes using X1 (one-versus-the-rest)
from sklearn.svm import LinearSVC
svm = LinearSVC(C=100,max_iter=50000)
svm.fit(X1, labels1)

# General labels using hyperplanes and X2
labels2 = svm.predict(X2)
print(np.unique(labels2))

# Perform differential expression
np.random.seed(0)
results = {}
start = time.time()

''' One-versus-one
for i, (c1, c2) in enumerate(itertools.combinations(np.unique(labels2), 2)):
    if (c1, c2) in results:
        continue
    #p_t = ttest_ind(X1[labels1 == c1].todense(), X1[labels1 == c2].todense())[1]
    p_t = ttest_ind(X1[labels1 == c1], X1[labels1 == c2])[1]
    p_t[np.isnan(p_t)] = 1
    #y = np.array(X2[labels2 == c1].todense())
    #z = np.array(X2[labels2 == c2].todense())
    y = np.array(X2[labels2 == c1])
    z = np.array(X2[labels2 == c2])
    #a = np.array(svm.coef_[i].todense()).reshape(-1)
    a = np.array(svm.coef_[i]).reshape(-1)
    b = svm.intercept_[i]
    p_tn, likelihood = tn.tn_test(y, z, a=a, b=b,
                                  learning_rate=1.,
                                  eps=1e-2,
                                  verbose=True,
                                  return_likelihood=True,
                                  num_iters=100000,
                                  num_cores=12)
    results[(c1, c2)] = (p_t, p_tn)
    print('c1: %5s\tc2: %5s\ttime elapsed: %.2fs'%(c1, c2, time.time()-start))
'''

''' One versus-the-rest
'''
for i, c1 in enumerate(np.unique(labels2)):
    #p_t = ttest_ind(X1[labels1 == c1].todense(), X1[labels1 == c2].todense())[1]
    p_t = ttest_ind(X1[labels1 == c1], X1[labels1 != c1])[1]
    p_t[np.isnan(p_t)] = 1
    #y = np.array(X2[labels2 == c1].todense())
    #z = np.array(X2[labels2 == c2].todense())
    y = np.array(X2[labels2 == c1])
    z = np.array(X2[labels2 != c1])
    #a = np.array(svm.coef_[i].todense()).reshape(-1)
    a = np.array(svm.coef_[i]).reshape(-1)
    b = svm.intercept_[i]
    p_tn, likelihood = tn.tn_test(y, z, a=a, b=b,
                                  learning_rate=1.,
                                  eps=1e-2,
                                  verbose=True,
                                  return_likelihood=True,
                                  num_iters=100000,
                                  num_cores=12)
    results[(c1)] = (p_t, p_tn)
    print('c1: %5s\ttime elapsed: %.2fs'%(c1, time.time()-start))

# open a file, where you ant to store the data
file = open('./results/ifnb_TN_test.pickle', 'wb')

# dump information to that file
pickle.dump(results, file)

# close the file
file.close()

# open a file, where you stored the pickled data
#file = open('./results/ifnb_TN_test.pickle', 'rb')

# dump information to that file
#results = pickle.load(file)

# close the file
#file.close()

# Result tidy
for i, (c1, c2) in enumerate(itertools.combinations(range(0,14), 2)):
    results[(c1,c2)] = results[(c1,c2)][1]
    
for i, c1 in enumerate(range(0,14)):
    results[(c1)] = results[(c1)][1]

p_frame = pd.DataFrame(results)
p_frame.to_csv("./results/ifnb_TN_test.csv",index = False,header=False)