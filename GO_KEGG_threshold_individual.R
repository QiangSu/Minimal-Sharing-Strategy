library(clusterProfiler)
library(org.Mm.eg.db)
library(DOSE)

# Step 1: Read the CSV file
file_path <- "C:/Users/Qiang/Desktop/SIAT/data/DEseq2_musHD_musWT_STR_result.csv"
degs <- read.csv(file_path, header = TRUE, sep = ",")

# Extract base name from the input file path
input_base_name <- tools::file_path_sans_ext(basename(file_path))

# Verify the data is read correctly
print(head(degs))

# Step 2: Check and print column names
print(colnames(degs))

# Step 3: Identify the correct columns for gene_id, log2FoldChange, and pvalue
gene_id_col <- "tracking_id"  # Adjust this column name to match your data
log2fc_col <- "log2FoldChange"  # Adjust this column name to match your data
pvalue_col <- "pvalue"  # Adjust this column name to match your data

if(!(gene_id_col %in% colnames(degs) & log2fc_col %in% colnames(degs) & pvalue_col %in% colnames(degs))) {
  stop("The necessary columns (tracking_id, log2FoldChange, pvalue) are not present in the data frame.")
}

# Step 4: Filter data based on p-value and log2FoldChange thresholds
filtered_degs <- degs[degs[[pvalue_col]] < 0.05 & (degs[[log2fc_col]] > -1 | degs[[log2fc_col]] < 1), ]

# Ensure filtered data is not empty
if (nrow(filtered_degs) == 0) {
  stop("No genes passed the filtering criteria (pvalue < 0.05, log2FoldChange > 1 or < -1). Please adjust the thresholds.")
}

# Step 5: Create gene list for enrichment analysis
gene_list <- setNames(filtered_degs[[log2fc_col]], filtered_degs[[gene_id_col]])
print(head(gene_list))

# Step 6: Convert gene identifiers to ENTREZID
gene_df <- bitr(names(gene_list), fromType = "ENSEMBLTRANS", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# Ensure that the conversion was successful
if (nrow(gene_df) == 0) {
  stop("Gene ID conversion failed. Please check your gene identifiers.")
}

# Step 7: Update gene list to use ENTREZID
gene_list <- gene_list[gene_df$ENSEMBLTRANS]
names(gene_list) <- gene_df$ENTREZID

# Step 8: Perform GO Enrichment Analysis for BP, MF, and CC separately
go_categories <- c("BP", "MF", "CC")
plot_list <- list()  # Create an empty list to store plots

for (go_cat in go_categories) {
  ego <- enrichGO(gene = names(gene_list),
                  OrgDb = org.Mm.eg.db,
                  keyType = "ENTREZID",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  ont = go_cat)
  
  # Check if any GO terms were enriched
  if (is.null(ego) || length(ego) == 0) {
    cat(paste("No significantly enriched GO terms found for", go_cat, ". Try less stringent thresholds or check your gene identifiers.\n"))
  } else {
    # View Results
    print(summary(ego))
    
    # Generate and store the barplot and dotplot
    barplot_ego <- barplot(ego, showCategory=40, title=paste("Top GO terms in", go_cat))
    dotplot_ego <- dotplot(ego, showCategory=40, title=paste("Top GO terms in", go_cat))
    
    # Store plots in the list for later display
    plot_list[[paste(go_cat, "barplot", sep = "_")]] <- barplot_ego
    plot_list[[paste(go_cat, "dotplot", sep = "_")]] <- dotplot_ego
    
    # Optionally, you can display the plots here in R Studio or Jupyter
    print(barplot_ego)
    print(dotplot_ego)
    
    # Save the results for each GO category, attaching the input file base name dynamically
    output_file_go <- paste0("C:/Users/Qiang/Desktop/SIAT/data/GO_enrichment_results_", go_cat, "_", input_base_name, ".csv")
    write.csv(as.data.frame(ego), file = output_file_go)
  }
}

# Display all six diagrams for BP, MF, and CC
par(mfrow=c(3, 2))  # Set up the plotting area for 6 plots (3 rows, 2 columns)

for (i in names(plot_list)) {
  plot(plot_list[[i]])  # Display each stored plot
}

# Step 9: Perform KEGG Enrichment Analysis
ekegg <- enrichKEGG(gene = names(gene_list),
                    organism = 'mmu',
                    keyType = "kegg",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2)

# Check if any KEGG pathways were enriched
if (is.null(ekegg) || nrow(as.data.frame(ekegg)) == 0) {
  cat("No significantly enriched KEGG pathways found. Try less stringent thresholds or check your gene identifiers.\n")
} else {
  # View Results
  print(as.data.frame(ekegg))  # Replace summary(ekegg) with as.data.frame(ekegg)
  
  # Generate barplot and dotplot for KEGG pathways
  barplot(ekegg, showCategory=10, title="Top KEGG pathways")
  dotplot(ekegg, showCategory=10, title="Top KEGG pathways")
  
  # Save the Results, attaching the input file base name dynamically
  output_file_kegg <- paste0("C:/Users/Qiang/Desktop/SIAT/data/KEGG_enrichment_results_", input_base_name, ".csv")
  write.csv(as.data.frame(ekegg), file = output_file_kegg)
}

