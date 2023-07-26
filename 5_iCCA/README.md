# Analysis of iCCA samples from 14 patients

- 5.1 Pre-clustering. Count matrix should be first downloaded from xxx and placed in this folder.
- 5.2 Run Festem. Cost 27 minutes (wall clock time) for each batch with 20 cores.
- 5.3 Run competing feature selection methods.
- 5.4 Cluster cells based on genes derived from different feature selection methods.
- 5.5 Based on clustering in 5.4, assigns DEGs to the clusters as their markers by the Scott-Knott test.
- 5.6 Figure 5, Figure S15-S19.


## Reference

Song, G., Shi, Y., Meng, L. et al. Single-cell transcriptomic analysis suggests two molecularly distinct subtypes of intrahepatic cholangiocarcinoma. Nat Commun 13, 1642 (2022). https://doi.org/10.1038/s41467-022-29164-0

## Prerequisite

- Download the preprocessed data from xxx

## Dependency

Apart from packages listed in "1_Simulation", the following R packages are also needed:

SeuratData==0.2.1, plyr==1.8.6, pROC==1.17.0.1, ScottKnott==1.3-0.