library(Seurat)
library(SeuratData)
library(ggplot2)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if a file path has been provided
if (length(args) < 3) {
  stop("Please provide a path to CSV file, metadata file, housekeeping gene file, and output path")
}

# Assign the file paths to variables
file_path <- args[1]
metadata_file_path <- args[2]
output_path <- args[3]

# Load count data
data_counts <- read.csv(file_path, row.names=1, header=T)
# Remove genes with zero counts across all cells
data_counts <- data_counts[rowSums(data_counts) > 0, ]

seurat_obj <- CreateSeuratObject(counts=data_counts)

# Load metadata
metadata <- read.csv(metadata_file_path, row.names = 1)
seurat_obj <- AddMetaData(seurat_obj, metadata = metadata)

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)

# Create pseudobulk for each cell type
pseudobulk <- AverageExpression(seurat_obj, 
    assays = "RNA", 
    group.by = "cell_type", 
    return.seurat = TRUE)
pseudobulk <- FindVariableFeatures(pseudobulk, selection.method = "vst", nfeatures = 2000)
pseudobulk <- RunPCA(pseudobulk, npcs = 8)
pseudobulk <- FindNeighbors(pseudobulk, dims = 1:8)
pseudobulk <- FindClusters(pseudobulk, resolution = 0.5)


expression_data <- GetAssayData(pseudobulk, assay = "RNA", layer = "data")
top_genes <- VariableFeatures(pseudobulk)[1:100]
expression_subset <- expression_data[top_genes, ]
expression_subset <- expression_subset[rowSums(expression_subset) > 0, ]

# K-means clustering
set.seed(42)
k <- 4
kmeans_result <- kmeans(t(expression_subset), centers = k)
pseudobulk$kmeans_cluster <- as.factor(kmeans_result$cluster)
heatmap <- DoHeatmap(pseudobulk, 
    features = top_genes, 
    group.by = "kmeans_cluster", 
    label = TRUE, 
    draw.lines = FALSE) + 
    ggtitle("K-Means Clustering of Pseudobulk Samples") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(axis.text.y = element_blank(),    # Optionally remove gene labels
    legend.position = "right",          # Optionally remove legend
    panel.grid = element_blank())  # Remove grid lines

# Save heatmap plot
ggsave(filename = file.path(output_path, "pseudobulk_kmeans_clustering_heatmap4.png"), 
       width = 10, height = 8)
