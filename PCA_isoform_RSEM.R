# Load necessary libraries
library(dplyr)
library(DESeq2)
library(ggplot2)
library(stringr)

# Change the path to the actual path of your CSV file
counts_path <- "C:/Users/Qiang/Desktop/SIAT/data/musWT_HD_RSEM_nobiascorrect_combined_est_counts.csv"

# Read in the data from the CSV file
count_data <- read.csv(counts_path, header = TRUE, sep = ",", check.names = FALSE)

# Ensure the selected columns exist
required_columns <- c("transcript_id", "musHD_STR-1_expected_count", "musHD_STR-2_expected_count", 
                      "musHD_STR-3_expected_count", "musWT_STR-1_expected_count", 
                      "musWT_STR-2_expected_count", "musWT_STR-3_expected_count")
missing_columns <- setdiff(required_columns, names(count_data))
if (length(missing_columns) > 0) {
  stop(paste("Missing required columns:", paste(missing_columns, collapse = ", ")))
}

# Select relevant columns for analysis
selected_count <- count_data %>%
  select(transcript_id,
         `musHD_STR-1_expected_count`,
         `musHD_STR-2_expected_count`,
         `musHD_STR-3_expected_count`,
         `musWT_STR-1_expected_count`,
         `musWT_STR-2_expected_count`,
         `musWT_STR-3_expected_count`)

# Filtering the counts (excluding the first column: transcript_id)
selected_count <- selected_count[rowMeans(selected_count[-1]) > 1,]

# Create count matrix and ensure no NA values
count_matrix <- as.matrix(round(selected_count[,-1]))  # Exclude the non-numeric 'transcript_id' column
rownames(count_matrix) <- selected_count$transcript_id  # Set 'transcript_id' as the row names

# Remove rows with NA values
count_matrix <- na.omit(count_matrix)

# Create sample information
data <- data.frame(
  Sample = c("musHD_STR-1_expected_count", "musHD_STR-2_expected_count", "musHD_STR-3_expected_count",
             "musWT_STR-1_expected_count", "musWT_STR-2_expected_count", "musWT_STR-3_expected_count"),
  Type = c("musHD", "musHD", "musHD", "musWT", "musWT", "musWT")
)
rownames(data) <- data$Sample
data$Sample <- NULL
data$Type <- as.factor(data$Type)

# Check if the column names of 'count_matrix' match the row names of 'data'
if (!all(rownames(data) %in% colnames(count_matrix))) {
  mismatched_names <- setdiff(rownames(data), colnames(count_matrix))
  stop(paste("Mismatch between sample names and count data columns:", paste(mismatched_names, collapse = ", ")))
}

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = data, design = ~ Type)

# Perform variance stabilizing transformation (vst) for PCA
vsd <- vst(dds, blind = FALSE)

# Plot PCA
pcaData <- plotPCA(vsd, intgroup = "Type", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = Type)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  coord_fixed() +  # Ensures equal scaling on both axes
  theme(
    axis.title.x = element_text(face = "bold", size = 10),  # Bold x-axis label with uniform size
    axis.title.y = element_text(face = "bold", size = 10),  # Bold y-axis label with uniform size
    axis.text.x = element_text(face = "bold", size = 10),   # Bold x-axis tick labels
    axis.text.y = element_text(face = "bold", size = 10),   # Bold y-axis tick labels
    plot.title = element_text(face = "bold", size = 10)     # Bold plot title with uniform size
  ) +
  ggtitle("PCA of musHD and musWT from RSEM")

# Output the PCA coordinates (PC1 and PC2) to a text file
output_file <- "C:/Users/Qiang/Desktop/SIAT/data/pca_coordinates_RSEM.txt"
write.table(pcaData, file = output_file, row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

# Confirm the file writing operation
cat("PCA coordinates have been written to:", output_file, "\n")

