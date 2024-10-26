library(Seurat)
library(ggplot2)
library(RColorBrewer) 

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if a file path has been provided
if (length(args) < 1) {
  stop("Please provide a path to CSV file and output path")
}

# Assign the file path to a variable
file_path <- args[1]
metadata_file_path <- args[2]
output_path <- args[3]

data_counts <- read.csv(file_path, row.names=1, header=T)
seurat_obj <- CreateSeuratObject(counts=data_counts)

# Preprocess
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)

# Load metadata
metadata <- read.csv(metadata_file_path, row.names = 1)
seurat_obj <- AddMetaData(seurat_obj, metadata = metadata)

# Run PCA
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))

# Run t-SNE
seurat_obj <- RunTSNE(seurat_obj, dims = 1:10)

# Plot t-SNE
custom_colors <- brewer.pal(n = 9, name = "Set1")
png(file.path(output_path, "tSNE_all_cells1.png"), width = 1000, height = 800)
DimPlot(seurat_obj, reduction = "tsne", group.by = "cell_type",pt.size = 4, cols = custom_colors)  
dev.off()



