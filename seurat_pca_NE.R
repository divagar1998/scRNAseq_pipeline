library(Seurat)
library(ggplot2)
library(dplyr)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if a file path has been provided
if (length(args) < 4) {
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

# Calculate coefficient of variation (CV) for each gene
cv_data <- apply(data_counts, 1, function(x) {
    mean_tpm <- mean(x[x > 0])  
    sd_tpm <- sd(x[x > 0])    
    cv <- sd_tpm / (mean_tpm + 1)  # Avoid division by zero
    return(cv)
})

cv_df <- data.frame(gene = rownames(data_counts), cv = cv_data)
# Remove NA values from cv_df
cv_df <- cv_df[!is.na(cv_df$cv), ]

# Calculate mean expression for each gene
cv_df$mean_expression <- rowMeans(data_counts[cv_df$gene, ]) + 1  # Adding 1 to avoid log(0)

# Plot mean expression vs CV
png(file.path(output_path, "mean_expression_vs_cv.png"), width = 1000, height = 800)
ggplot(cv_df, aes(x = mean_expression, y = cv)) + geom_point(alpha = 0.5) +
    scale_x_log10() + 
    labs(title = "Mean Expression vs Coefficient of Variation (CV)", x = "Mean Expression (log scale)", y = "Coefficient of Variation (CV)")
dev.off()

# Create polynomial model for mean expression vs CV
polynomial_fit <- lm(cv ~ poly(log(mean_expression), 2, raw = TRUE), data = cv_df)

# Generate predictions for plotting
cv_df$predicted_cv <- predict(polynomial_fit, newdata = cv_df)
cv_df$predicted_cv[cv_df$predicted_cv < 0] <- 0.001

# Plot mean expression vs CV with fitted polynomial curve
png(file.path(output_path, "mean_expression_vs_cv_with_fit4.png"), width = 1000, height = 800)
ggplot(cv_df, aes(x = mean_expression, y = cv)) +
    geom_point(alpha = 0.5) +
    geom_line(aes(y = predicted_cv), color = "blue", linewidth = 1) + 
    scale_x_log10() +  
    labs(title = "Mean Expression vs Coefficient of Variation (CV) with Polynomial Fit", 
       x = "Mean Expression (log scale)", 
       y = "Coefficient of Variation (CV)")
dev.off()

# Calculate the ratio of measured CV to predicted CV
cv_df$cv_ratio <- cv_df$cv / cv_df$predicted_cv
cv_df$cv_ratio[is.na(cv_df$cv_ratio)] <- NA 

# Plot distribution of CV ratio
png(file.path(output_path, "cv_ratio_distribution2.png"), width = 1000, height = 800)
ggplot(cv_df, aes(x = log10(cv_ratio + 1))) +  # Add 1 to avoid log(0)
    geom_histogram(binwidth = 0.01, fill = "lightblue", color = "black", alpha = 0.7) +
    labs(title = "Distribution of CV Ratio (Log10 Scale)", x = "Log10(CV Ratio + 1)", y = "Frequency") +
    xlim(c(0, max(log10(cv_df$cv_ratio + 1), na.rm = TRUE))) # Limit x-axis if needed
dev.off()

# Take the top 500 genes with highest CV ratio
top_500_genes <- cv_df[order(cv_df$cv_ratio, decreasing = TRUE), ][1:500, ]

# Get the top 100 genes with the highest expression and CV ratio > 1.5
# Calculate mean expression for each gene
cv_df$mean_expression <- rowMeans(data_counts[rownames(cv_df), ], na.rm = TRUE)
# Filter genes based on CV ratio and mean expression
top_100_genes <- cv_df[cv_df$cv_ratio > 1.5, ]
top_100_genes <- top_100_genes[order(top_100_genes$mean_expression, decreasing = TRUE), ][1:100, ]

# Combine the selected genes
selected_genes <- unique(c(top_500_genes$gene, top_100_genes$gene))
# Check if selected genes are in the data
valid_genes <- selected_genes[selected_genes %in% rownames(data_counts)]

# Preprocess
seurat_obj <- NormalizeData(seurat_obj)
# Subset the Seurat object to include only the selected genes
seurat_selected <- CreateSeuratObject(counts = data_counts[valid_genes, ])
seurat_selected <- AddMetaData(seurat_selected, metadata = metadata)
seurat_selected <- NormalizeData(seurat_selected)
seurat_selected <- FindVariableFeatures(seurat_selected)
seurat_selected <- ScaleData(seurat_selected)

# Run PCA
seurat_selected <- RunPCA(seurat_selected, features = VariableFeatures(seurat_selected))

# Plot PCA
png(file.path(output_path, "pca_selected_genes.png"), width = 1000, height = 800)
DimPlot(seurat_selected, reduction = "pca", group.by = "cell_type", pt.size = 4) + 
    labs(title = "PCA of Selected Genes", x = "PC1", y = "PC2")
dev.off()

# Plot PCA with gene expression colour coding
gene_of_interest <- "Chga:NM-007693.1"
gene_expression <- FetchData(seurat_obj, vars = gene_of_interest)
gene_expression_vector <- as.vector(gene_expression[, 1]) 

# Assign the gene expression values to seurat_selected for the selected cells
cell_indices <- match(colnames(seurat_selected), colnames(seurat_obj))
seurat_selected$gene_expression <- gene_expression_vector[cell_indices]

# Specify the cell types to highlight
cell_types_to_highlight <- c("NE_AT2", "NE_ciliated", "NE_club", "NE_stromal")  
# Create a new column in seurat_selected to indicate whether a cell is of the specified types
seurat_selected$highlight <- ifelse(seurat_selected$cell_type %in% cell_types_to_highlight, "Highlight", "Normal")

# Extract PCA coordinates from the main Seurat object
pca_coords <- Embeddings(seurat_selected, "pca")
# Add PCA coordinates to the selected object for plotting
seurat_selected$pca_1 <- pca_coords[, 1]
seurat_selected$pca_2 <- pca_coords[, 2]

png(file.path(output_path, "feature_plot_chga_circled_pca_1.png"), width = 1000, height = 800)
FeaturePlot(seurat_selected, features = "gene_expression", cols = c("lightgrey", "red"), pt.size = 4) + 
    geom_point(data = subset(seurat_selected@meta.data, highlight == "Highlight"), 
    aes(x = pca_1, y = pca_2), 
    shape = 1, size = 5, color = "black", stroke = 1) + 
    labs(title = paste("Feature Plot of", gene_of_interest, "with Transitioning NE cells Circled"))
dev.off()




