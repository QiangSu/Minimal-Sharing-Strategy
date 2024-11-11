# Load necessary libraries
library(dplyr)
library(DESeq2)
library(stringr)
library(ggplot2)  # For PCA plotting

# Change the path to the actual path of your CSV file
counts_path <- "C:/Users/Qiang/Desktop/SIAT/data/combined_results_GaussF_GC_mus_count.csv"

# Read in the data from the CSV file
count_data <- read.csv(counts_path, header = TRUE, sep = ",", check.names = FALSE)

# Preprocess "Present_in_Transcripts" to ensure consistency by trimming whitespace and converting to uppercase
count_data <- count_data %>%
  mutate(Present_in_Transcripts = as.character(Present_in_Transcripts),  # Convert to character to ensure proper manipulation
         Present_in_Transcripts = str_trim(Present_in_Transcripts),      # Remove leading/trailing whitespace
         Present_in_Transcripts = toupper(Present_in_Transcripts))       # Convert to uppercase for uniformity

# Ensure there are unique rows based on "Present_in_Transcripts"
count_data <- count_data %>%
  distinct(Present_in_Transcripts, .keep_all = TRUE)

# Make sure that the selected columns exist
required_columns <- c("Present_in_Transcripts", 
                      "results_file_GC_musHD_STR-1", "results_file_GC_musHD_STR-2", "results_file_GC_musHD_STR-3", 
                      "results_file_GC_musWT_STR-1", "results_file_GC_musWT_STR-2", "results_file_GC_musWT_STR-3")
missing_columns <- setdiff(required_columns, names(count_data))
if (length(missing_columns) > 0) {
  stop(paste("Missing required columns:", paste(missing_columns, collapse = ", ")))
}

# Select relevant columns for analysis
selected_count <- count_data %>%
  select(Present_in_Transcripts,
         `results_file_GC_musHD_STR-1`,
         `results_file_GC_musHD_STR-2`,
         `results_file_GC_musHD_STR-3`,
         `results_file_GC_musWT_STR-1`,
         `results_file_GC_musWT_STR-2`,
         `results_file_GC_musWT_STR-3`)

# Filtering the counts based on row means (excluding the 'Present_in_Transcripts' column)
selected_count <- selected_count[rowMeans(selected_count[-1]) > 1,]

# Additional filtering to keep genes expressed in at least 5 samples (since 6 samples total)
#selected_count <- selected_count[rowSums(selected_count[-1] > 0) >= 5,]

# Convert the count matrix to a proper format
count_matrix <- as.matrix(round(selected_count[,-1]))  # Exclude 'Present_in_Transcripts'
rownames(count_matrix) <- selected_count$Present_in_Transcripts  # Set gene names as row names

# Remove rows with NA values
count_matrix <- na.omit(count_matrix)

# Create the sample information data frame directly in the script
data <- data.frame(
  Sample = c("results_file_GC_musHD_STR-1", "results_file_GC_musHD_STR-2", "results_file_GC_musHD_STR-3", 
             "results_file_GC_musWT_STR-1", "results_file_GC_musWT_STR-2", "results_file_GC_musWT_STR-3"),
  Type = c("musHD", "musHD", "musHD", "musWT", "musWT", "musWT")
)
rownames(data) <- data$Sample
data$Sample <- NULL

# Set 'Type' as a factor
data$Type <- as.factor(data$Type)

# Create the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = data,
                              design = ~ Type)

# Perform variance stabilizing transformation (VST) for PCA
vsd <- vst(dds, blind = FALSE)

# PCA plot using DESeq2's plotPCA function
pcaData <- plotPCA(vsd, intgroup = "Type", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Updated PCA plot using ggplot2 with requested settings
p <- ggplot(pcaData, aes(PC1, PC2, color = Type)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", round(percentVar[1], 2), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2], 2), "% variance")) +
  theme_minimal() +  # Apply minimal theme for cleaner design
  coord_fixed() +    # Ensures equal scaling on both axes (square plot)
  theme(
    axis.title.x = element_text(face = "bold", size = 10),  # Bold x-axis label with larger size
    axis.title.y = element_text(face = "bold", size = 10),  # Bold y-axis label with larger size
    axis.text.x = element_text(face = "bold", size = 10),   # Bold x-axis tick labels
    axis.text.y = element_text(face = "bold", size = 10),   # Bold y-axis tick labels
    plot.title = element_text(face = "bold", size = 10)     # Adjust plot title to size 10
  ) +
  ggtitle("PCA of musHD and musWT from GaussF")  # Title with your specific formatting

# Print the updated PCA plot
print(p)

# Output the PCA coordinates (PC1 and PC2) to a text file
output_file <- "C:/Users/Qiang/Desktop/SIAT/data/pca_coordinates_GaussF.txt"
write.table(pcaData, file = output_file, row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

# Confirm the file writing operation
cat("PCA coordinates have been written to:", output_file, "\n")

