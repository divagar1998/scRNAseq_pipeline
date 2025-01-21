# install these packages in env
# python/3.12.0 
# scanpy/1.10.4
# anndata/0.11.0
# pandas/2.2.3
# scvi/1.2.1
# scikit-misc/0.5.1

import scanpy as sc
import anndata as ad 
import pandas as pd 
import scvi

count_matrix_path = "/hpc/users/divagt01/watanabe/Divagar/transcriptomic_data/Tang/OMIX002441-01.csv"

# for scanpy the genes are columns and cells are rows, vice versa for Seurat
df = pd.read_csv(count_matrix_path, index_col=0, header=0)
df_transposed = df.T

# load dataframe into scanpy object
tang_adata = sc.AnnData(df_transposed.values)
tang_adata.obs_names = df_transposed.index
tang_adata.var_names = df_transposed.columns

# filter out genes that are not present in at least 100 cells
sc.pp.filter_genes(tang_adata, min_cells=100)
# only keep the top 2000 variable genes
sc.pp.highly_variable_genes(tang_adata, n_top_genes=2000, subset=True, flavor="seurat_v3")

# doublet removal should be done per sample before integration
# scvi.model.SCVI.setup_anndata(tang_adata)
# vae = scvi.model.SCVI(tang_adata)
# vae.train()

# checking for and removing mitochondrial genes
