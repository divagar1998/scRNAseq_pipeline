library(Seurat)
library(ggplot2)
library(slingshot)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the necessary file paths have been provided
if (length(args) < 3) {
  stop("Please provide paths to CSV file, metadata file and output path")
}

# Assign the file paths to variables
file_path <- args[1]
metadata_file_path <- args[2]
output_path <- args[3]

# Load count data
data_counts <- read.csv(file_path, row.names=1, header=T)

# Remove genes with zero counts across all cells
data_counts <- data_counts[rowSums(data_counts) > 0, ]

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts=data_counts)

# Load metadata
metadata <- read.csv(metadata_file_path, row.names = 1)
seurat_obj <- AddMetaData(seurat_obj, metadata = metadata)

# Normalize, identify variable features, scale, and run PCA
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)

# Clustering the cells
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj)

# Extract PCA coordinates and clusters
pca_coords <- Embeddings(seurat_obj, "pca")
clusters <- seurat_obj$seurat_clusters

# Run Slingshot
slingshot_obj <- slingshot(pca_coords, clusterLabels = clusters, reducedDim = "PCA")

# Extract pseudotime lineage data
lineages <- slingPseudotime(slingshot_obj)

# Prepare data for ggplot
trajectory_df <- as.data.frame(pca_coords)
trajectory_df$cluster <- as.factor(clusters)

# Convert lineages to data frame and ensure it has appropriate column names
lineage_df <- data.frame(PC_1 = pca_coords[lineages[, 1], 1], PC_2 = pca_coords[lineages[, 1], 2],
                          pseudotime = lineages[, 1])
                          colnames(lineage_df) <- paste0("PC_", seq_len(ncol(lineage_df)))
# Plot with ggplot
p <- ggplot(trajectory_df, aes(x = PC_1, y = PC_2, color = cluster)) +
    geom_point() +
    geom_line(data = lineage_df, aes(x = PC_1, y = PC_2), color = "black", linewidth = 1) +
    labs(title = "Slingshot Trajectory Plot")

# Save the plot to a file
ggsave(filename = file.path(output_path, "trajectory_plot.png"), plot = p, width = 8, height = 6)
