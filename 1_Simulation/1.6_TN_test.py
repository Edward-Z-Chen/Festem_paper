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
import tracemalloc

os.chdir("./results")
time_cost = np.zeros(20)
tn_result = list()
memory_usage = list()

# Louvain
def scanpy_cluster(X, features, tsne=None, plot=False, resolution=1.0):
    "Scanpy implementation of Seurat's pipeline"

    adata = sc.AnnData(X=X)
    adata.var['genes_ids'] = features

    # preprocessing
    # sc.pp.filter_cells(adata, min_genes=200)
    # sc.pp.filter_genes(adata, min_cells=3)
    # mito_genes = adata.var_names.str.startswith('MT-')
    # adata.obs['percent_mito'] = np.sum(
    #    adata[:, mito_genes].X, axis=1) / np.sum(adata.X, axis=1)
    # adata.obs['n_counts'] = adata.X.sum(axis=1)

    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=1000, flavor="seurat_v3")
    adata = adata[:, adata.var['highly_variable']]
    # sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])
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

for batch in range(1,21):
    # Load data
    adata = read_r("NB_200DE_2type_"+str(batch)+".RData")
    #adata = read_r("NB_400DE_2type_"+str(batch)+".RData")
    #adata = read_r("NB_200DE_5type_"+str(batch)+".RData")
    #adata = read_r("NB_400DE_5type_"+str(batch)+".RData")
    
    counts = adata["counts"]
    features = counts.index.values
    counts = counts.transpose()
    counts = counts.to_numpy()
    counts = counts.astype(int)
    del adata

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
    for k in range(1,101):
        labels1 = scanpy_cluster(X1, features, plot=False, resolution=0.01 * k)
        if len(np.unique(labels1)) == 2:
            break
        # If five groups    
        #if len(np.unique(labels1)) == 5:
        #    break

    # Fit hyperplanes using X1
    from sklearn.svm import SVC
    svm = SVC(kernel='linear', C=100)
    svm.fit(X1, labels1)

    # General labels using hyperplanes and X2
    labels2 = svm.predict(X2)
    print("Number of clusters: {}".format(len(np.unique(labels2))))
    # Perform differential expression
    np.random.seed(4321)
    results = {}
    start = time.time()
    tracemalloc.start()
    
    if len(np.unique(labels2)) == 2:
        y = np.array(X2[labels2 == 0])
        z = np.array(X2[labels2 != 0])
        a = np.array(svm.coef_[0]).reshape(-1)
        b = svm.intercept_[0]
        p_tn, likelihood = tn.tn_test(y, z, a=a, b=b,
                                      learning_rate=1.,
                                      eps=1e-2,
                                      verbose=True,
                                      return_likelihood=True,
                                      num_iters=100000,
                                      num_cores=12)
        results[(0)] = p_tn
        print('c1: %5s\ttime elapsed: %.2fs'%(0, time.time()-start))
    else:
        for i, c1 in enumerate(np.unique(labels2)):
            #p_t = ttest_ind(X1[labels1 == c1].todense(), X1[labels1 == c2].todense())[1]
            #p_t = ttest_ind(X1[labels1 == c1], X1[labels1 != c1])[1]
            #p_t[np.isnan(p_t)] = 1
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
            results[(c1)] = p_tn
            print('c1: %5s\ttime elapsed: %.2fs'%(c1, time.time()-start))
        
    tmp = tracemalloc.get_traced_memory()
    memory_usage.append(tmp)
    tracemalloc.stop()
    time_cost[batch-1] = (time.time() - start)*12
    # for i, c1 in enumerate(np.unique(labels2)):
    #     results[(c1)] = results[(c1)][1]
    p_frame = pd.DataFrame(results)
    tn_result.append(p_frame)
    
    
#p_frame.to_csv("pbmc_TN_test_one_versus_the_rest.csv",index = False,header=False)

# For other settings, change the following file names
np.save("NB_200DE_2celltype_TN_test_time.npy",time_cost)
# open a file, where you ant to store the data
file = open('NB_200DE_2celltype_TN_test', 'wb')

# dump information to that file
pickle.dump(tn_result, file)

# close the file
file.close()

file = open('NB_200DE_2celltype_TN_test_memory', 'wb')

# dump information to that file
pickle.dump(memory_usage, file)

# close the file
file.close()

# Load in results and write csv
file = open('NB_200DE_2celltype_TN_test', 'rb')
results = pickle.load(file)
file.close()
p_frame = np.zeros((20000,20))
for j in range(0,20):
    #p_frame[:,j] = results[j].iloc[:,0]
    p_frame[:,j] = np.min(results[j],axis = 1) * results[j].shape[1]
p_frame = pd.DataFrame(p_frame)
p_frame.to_csv("NB_200DE_2celltype_TN_test.csv",index = False,header=False,na_rep = "NA")

# The second column is peak
file = open('NB_200DE_2celltype_TN_test_memory', 'rb')
memory_usage = pickle.load(file)
file.close()
memory_usage = pd.DataFrame(memory_usage)
memory_usage = memory_usage/(1024**2)
memory_usage.to_csv("NB_200DE_2celltype_TN_test_memory.csv",index = False,header=False,na_rep = "NA")

time_cost = np.load("NB_200DE_2celltype_TN_test_time.npy")
time_cost = pd.DataFrame(time_cost)
time_cost.to_csv("NB_200DE_2celltype_TN_test_time.csv",index = False,header=False,na_rep = "NA")
