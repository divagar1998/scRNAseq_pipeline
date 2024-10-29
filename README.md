# scRNAseq_pipeline
This repository consists of code used to analyze the scRNA-seq data generated from Ouadah,...,Krasnow doi: 10.1016/j.cell.2019.09.010. and Ireland,...,Oliver doi: 10.1016/j.ccell.2020.05.001. \

Ouadah et. al has scRNA-seq data from mouse Pulmonary Neuroendocrine Cells (PNECs) that are reprogramming into non-NE cell types upon napthelene injury \
Ireland et. al has scRNA-seq data from mouse models of SCLC that are transitioning from NE-high to NE-low SCLC \

## SEURAT
`seurat_tsne.R` uses Seurat to generate a tSNE plot of all the cell types in mouse lung after napthelene injury in Ouadah et. al \
`seurat_pca_NE.R` uses Seurat to generate PCA plots of NE cells in mouse lung after napthelene injury in Ouadah et. al. Feature plot allows me to show expression of certain genes in the PCA. \
`convert_count_to_tpm.R` converts kallisto pseudo counts to TPM \


