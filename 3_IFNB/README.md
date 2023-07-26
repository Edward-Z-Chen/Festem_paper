# Analysis of the IFNB dataset (Kang dataset)

- 3.1 Download preprocessed data through *SeuratData* and pre-clustering for the control group. Unless otherwise stated, the analyses are for the control group. The downloaded data is preprocessed and annotated according to https://satijalab.org/seurat/v3.0/immune_alignment.html.
- 3.2 Run Festem. Cost 69 seconds (wall clock time) with 12 cores.
- 3.3 Run competing DEG detection methods.
- 3.4 Run TN test.
- 3.5 Run competing feature selection methods.
- 3.6 Clustering and UMAP reduction based on top 2500 genes selected by each method and the top 25 PCs.
- 3.7 Calculate CH indices of the PCs and clustering results derived from genes selected by feature selection methods.
- 3.8 Construct silver standard of DEGs based on housekeeping genes downloaded from https://housekeeping.unicamp.br/Housekeeping_GenesHuman.RData.
- 3.9 Festem algorithm and clustering of the stimulated group with the same procedure.
- 3.10 Figure 3A (top right), Figure 4E-G, Figure S5 (second from top), Figure S7 (top right), Figure S10, Figure S11, Figure S12, Figure S13 and Figure S14. In Figure S13A and B, we performed integrated analysis of the control and stimulated groups with *Harmony*. P-values of features were combined with the Bonferroni procedure and FDR were controlled with the Benjamini-Hochberg procedure.

**Remark**: For results in 3.3 and 3.4, we only provide a tidy version named "ifnb_ctrl_DEG_results.RData" derived from 3.10 in the "results" folder.

## Reference

Kang, H.M. et al. Multiplexed droplet single-cell RNA-sequencing using natural genetic variation. Nature biotechnology 36, 89-94 (2018).

Vovk, V. & Wang, R. Combining P-Values Via Averaging. Political Methods: Quantitative Methods eJournal (2012).

## Prerequisite

- Download housekeeping gene list from https://housekeeping.unicamp.br/Housekeeping_GenesHuman.RData and place it in **this** folder.

## Dependency

Apart from packages listed in "1_Simulation", the following R packages are also needed:

SeuratData==0.2.1, Harmony==0.1.0, plyr==1.8.6, pROC==1.17.0.1, ggsci==2.9, cowplot==1.1.1, ggalluvial==0.12.3