# Load necessary library
library(dplyr)

count_matrix <- read.csv("/hpc/users/divagt01/watanabe/Divagar/krasnow_scrnaseq/kallisto_counts_subset.csv", row.names = 1)
# Remove zsGreen 
count_matrix <- count_matrix[!rownames(count_matrix) %in% "zsGreen", ]

# Multiply counts by 66 to get TPM
count_matrix_scaled <- count_matrix * 66

# Calculate total counts per column
total_counts <- colSums(count_matrix_scaled)

# Calculate TPM
tpm_matrix <- sweep(count_matrix_scaled, 2, total_counts, FUN = "/") * 1e6

# Set TPM values less than 10 to 0
tpm_matrix[tpm_matrix < 10] <- 0

# Convert to data frame
tpm_df <- as.data.frame(tpm_matrix)

# Save the TPM matrix
write.csv(tpm_df, "/hpc/users/divagt01/watanabe/Divagar/krasnow_scrnaseq/kallisto_counts_subset_tpm.csv", row.names = TRUE)
