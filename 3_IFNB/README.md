# Analysis of the IFNB dataset

- 3.1 Download preprocessed data through *SeuratData* and pre-clustering. The downloaded data is preprocessed and annotated according to https://satijalab.org/seurat/v3.0/immune_alignment.html.

## Reference
Kang, H.M. et al. Multiplexed droplet single-cell RNA-sequencing using natural genetic variation. Nature biotechnology 36, 89-94 (2018).

## Dependency

Besides the packages listed in "1_Simulation", the following R packages are also needed:

SeuratData==0.2.1, plyr==1.8.6, pROC==1.17.0.1, ggsci==2.9, cowplot==1.1.1, ggalluvial==0.12.3