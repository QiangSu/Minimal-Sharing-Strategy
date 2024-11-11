# Load necessary libraries
library(readr)
library(VennDiagram)   # Limited to 5 sets, better to use UpSetR for more
library(UpSetR)        # For UpSet plots, better for more than 5 sets
library(ggVennDiagram) # To handle more than 5 sets in a Venn-like fashion
library(dplyr)
library(ggplot2)
# Load the data using the correct file path
data <- read_csv('C:/Users/Qiang/Desktop/SIAT/data/methods_filtered_isoform_combined.csv')

# Confirm data is loaded correctly
glimpse(data)

# Extract the sets from each column, ensuring no NAs
GaussF <- data$GaussF[!is.na(data$GaussF)]
Salmon <- data$Salmon[!is.na(data$Salmon)]
Kallisto <- data$Kallisto[!is.na(data$Kallisto)]
RSEM <- data$RSEM[!is.na(data$RSEM)]
Cufflinks <- data$Cufflinks[!is.na(data$Cufflinks)]
StringTie <- data$StringTie[!is.na(data$StringTie)]

# Create a list of sets for visualization and further analysis
gene_sets <- list(
  GaussF = GaussF,
  Salmon = Salmon,
  Kallisto = Kallisto,
  RSEM = RSEM,
  Cufflinks = Cufflinks,
  StringTie = StringTie
)

# Option 1: Using ggVennDiagram for a more than 5 sets Venn-like visualization
# Customize the plot: Set individual colors and apply transparency, remove counts
ggVennDiagram(gene_sets, label = "count") +
  scale_fill_gradient(low = "#FFFFFF", high = "#0073C2") +
  scale_color_manual(values = c("#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF")) +
  theme(legend.position = "none") +
  theme_void() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    text = element_text(size = 12, face = "bold") # Increase label size and make bold
  )


# Option 2: Using UpSetR for a scalable and clear intersection visualization

library(UpSetR)
library(grid)

# Create the binary matrix as before
all_genes <- unique(unlist(gene_sets))
gene_matrix <- data.frame(
  GaussF = as.integer(all_genes %in% GaussF),
  Salmon = as.integer(all_genes %in% Salmon),
  Kallisto = as.integer(all_genes %in% Kallisto),
  RSEM = as.integer(all_genes %in% RSEM),
  Cufflinks = as.integer(all_genes %in% Cufflinks),
  StringTie = as.integer(all_genes %in% StringTie)
)

# Set gene IDs as row names
rownames(gene_matrix) <- all_genes

# Define colors for each method
colors <- c("GaussF" = "red", 
            "Salmon" = "blue", 
            "Kallisto" = "green", 
            "RSEM" = "purple", 
            "Cufflinks" = "orange", 
            "StringTie" = "pink")

# Create the UpSet plot with custom colors, order by frequency, and adjust font sizes
upset_plot <- upset(gene_matrix, 
                    sets = c("GaussF", "Salmon", "Kallisto", "RSEM", "Cufflinks", "StringTie"),
                    order.by = "freq",
                    sets.bar.color = colors,
                    text.scale = c(2.5, 2.5, 2.5, 2.5, 2.5, 2.5))  # Adjust the numeric values to control the font size

# Add bold formatting to the labels
grid.text <- function(..., gp = gpar(fontface = "bold")) {
  grid::grid.text(...)
}

# Saving the plot as a TIFF file
tiff("C:/Users/Qiang/Desktop/SIAT/data/upset_plot.tiff", 
     width = 20, height = 10, units = "in", res = 1000) # Adjust width and height as needed

# Print and display the UpSet plot
print(upset_plot)

# Add bold font for set names and other grid elements
#grid.text(label = "Gene Expression Methods", x = 0.5, y = 0.95, gp = gpar(fontsize = 20, fontface = "bold"))

# Close the TIFF device
dev.off()
