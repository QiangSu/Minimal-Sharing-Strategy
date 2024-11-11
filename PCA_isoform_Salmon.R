# Load necessary libraries
library(dplyr)
library(DESeq2)
library(ggplot2)
library(stringr)

# Change the path to the actual path of your CSV file
counts_path <- "C:/Users/Qiang/Desktop/SIAT/data/musWT_HD_salmon_nobiascorrect_combined_quant.csv"

# Read in the data from the CSV file
count_data <- read.csv(counts_path, header = TRUE, sep = ",", check.names = FALSE)

# Inspect loaded data to confirm it is as expected
head(count_data)


print(nrow(count_data))
print(head(count_data))

# The input columns will be changed according to the different comparing groups
required_columns <- c("Name", "musHD_STR-1_NumReads", "musHD_STR-2_NumReads", "musHD_STR-3_NumReads", "musWT_STR-1_NumReads", "musWT_STR-2_NumReads", "musWT_STR-3_NumReads")
missing_columns <- setdiff(required_columns, names(count_data))
if(length(missing_columns) > 0) {
  stop(paste("Missing required columns:", paste(missing_columns, collapse=", ")))
}

# Select Present_in_Transcripts, results_file_GC_musHD_COX-1, results_file_GC_musHD_COX-2, results_file_GC_musHD_COX-3, results_file_GC_musHD_STR-1,results_file_GC_musHD_STR-2,results_file_GC_musHD_STR-3 for analysis
# The input columns will be changed according to the different comparing groups
selected_count <- dplyr::select(count_data, Name,
                                `musHD_STR-1_NumReads`,
                                `musHD_STR-2_NumReads`,
                                `musHD_STR-3_NumReads`,
                                `musWT_STR-1_NumReads`,
                                `musWT_STR-2_NumReads`,
                                `musWT_STR-3_NumReads`)

# Print the resulting data frame to check the selection
print(head(selected_count))

# Filtering the counts
# Exclude the first column (Present_in_Transcripts column) when calculating row means
selected_count <- selected_count[rowMeans(selected_count[-1]) > 1,]

# Print the head of the selected_count data frame after filtering
print(nrow(selected_count))

# Ensure there are no NA values in the count matrix
count_matrix <- as.matrix(round(selected_count[,-1])) # Exclude the non-numeric 'Present_in_Transcripts' column
rownames(count_matrix) <- selected_count$Name  # Assuming 'Present_in_Transcripts' are the gene/feature identifiers

# Remove rows with NA values
count_matrix <- na.omit(count_matrix)

# Create the sample information data frame directly in the script
data <- data.frame(
  Sample = c("musHD_STR-1_NumReads", "musHD_STR-2_NumReads", "musHD_STR-3_NumReads", "musWT_STR-1_NumReads", "musWT_STR-2_NumReads", "musWT_STR-3_NumReads"),
  Type = c("musHD", "musHD", "musHD", "musWT", "musWT", "musWT")
)
rownames(data) <- data$Sample
data$Sample <- NULL

# Preview data to confirm proper load
print(head(data))

# Set 'Type' as a factor
data$Type <- as.factor(data$Type)

# Check to ensure that the column names of 'count_matrix' match the row names of 'data'
if (!all(rownames(data) %in% colnames(count_matrix))) {
  mismatched_names <- setdiff(rownames(data), colnames(count_matrix))
  print(mismatched_names)
  stop("Mismatch between sample names and selected count data columns")
}

# Create the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = data,
                              design = ~ Type)

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
  ggtitle("PCA of musHD and musWT from Salmon")

# Output the PCA coordinates (PC1 and PC2) to a text file
output_file <- "C:/Users/Qiang/Desktop/SIAT/data/pca_coordinates_salmon.txt"
write.table(pcaData, file = output_file, row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

# Confirm the file writing operation
cat("PCA coordinates have been written to:", output_file, "\n")
