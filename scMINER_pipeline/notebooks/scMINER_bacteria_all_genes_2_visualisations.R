### present working directory: C:\Users\natha\OneDrive\Bio-informatica_25-26\International Internship\Internship-25-26\scMINER_pipeline
### scMINER pipeline - bacteria_all_genes_combined
### BF (biofilm) + P (planktonic) dataset
### 7884 cells x 2853 genes
# ==============================================================================
# STEP 1: LOAD THE COMBINED DATASET
# ==============================================================================
# Both datasets (BF and P) were already filtered, quality-controlled and
# normalised in Python (log1p-CPM values confirmed), and already concatenated
# into a single combined dataset → QC, filtration and normalisation are skipped.
# $X   → expression matrix (cells x genes): 7884 x 2853
# $obs → cell metadata (sample)
# $var → gene metadata (empty, but var_names contains the gene names)

library(scMINER)
library(Matrix)
library(ggplot2)
library(anndata)

adata_combined <- readInput_h5ad(h5ad_file = "./DATA/bacteria_all_genes_combined.h5ad")

# Add prefix based on sample column (bf_ or p_)
sample_labels <- adata_combined$obs$sample  # "bf" or "p"
cell_numbers  <- adata_combined$obs_names   # "68", "96", etc.
adata_combined$obs_names <- paste0(sample_labels, "_", cell_numbers)

# Verify
print(head(adata_combined$obs_names, 5))  # Should show bf_68, bf_96...

# ==============================================================================
# STEP 2: BUILD SPARSE EXPRESSION MATRIX
# ==============================================================================
# Transpose: cells x genes → genes x cells (2853 x 7884)
# FIX: use var_names instead of rownames(var) because var is empty

combined_mtx <- as(t(adata_combined$X), "dgCMatrix")
rownames(combined_mtx) <- adata_combined$var_names   # ← 2853 gene names (dnaA, gyrB...)
colnames(combined_mtx) <- adata_combined$obs_names   # ← 7884 cell names (bf_68, p_96...)

print(dim(combined_mtx))  # expected: 2853 genes x 7884 cells

# ==============================================================================
# STEP 3: CREATE SparseEset OBJECT
# ==============================================================================

cell_metadata <- as.data.frame(adata_combined$obs)
rownames(cell_metadata) <- adata_combined$obs_names  # ← ensure row names match

combined.eset <- createSparseEset(input_matrix = combined_mtx,
                                  cellData     = cell_metadata,
                                  projectID    = "BaSSSh_bacteria",
                                  addMetaData  = TRUE)

dim(combined.eset)          # expected: 2853 x 7884
head(pData(combined.eset))  # verify sample column is present

# ==============================================================================
# STEP 4: EXPORT AS MICA INPUT FILE
# ==============================================================================

#generateMICAinput(input_eset  = combined.eset,
                  #output_file = "C:/Users/natha/OneDrive/Bio-informatica_25-26/International Internship/Internship-25-26/scMINER_pipeline/MICA_INPUT/MICA_all_genes_Input.txt")

# ==============================================================================
# STEP 5: RUN MICA (TERMINAL — NOT IN R)
# ==============================================================================
# Heavy computation — run the following command in your terminal, not here:
#
#   mica ge -i "./MICA/MICA_all_genes_Input.txt" \
#           -o "./MICA/micaOutput" \
#           -nw 6 -nc 6

# ==============================================================================
# STEP 6: AUTOMATED DIRECTORY LOOP FOR DIRECT SILHOUETTE SCORE CALCULATION
# ==============================================================================



# Ensure the 'cluster' package is installed for standard validation metrics
if (!requireNamespace("cluster", quietly = TRUE)) {
  install.packages("cluster")
}
library(cluster)

# Your exact project and MDS directory paths
base_mds_dir    <- "C:/Users/natha/OneDrive/Bio-informatica_25-26/International Internship/Internship-25-26/scMINER_pipeline/MICA_OUTPUT/mica_mds_output_all_genes"
celltypist_file <- "C:/Users/natha/OneDrive/Bio-informatica_25-26/International Internship/Internship-25-26/scMINER_pipeline/CELLTYPIST/celltypist_predictions.txt"
results_dir     <- "C:/Users/natha/OneDrive/Bio-informatica_25-26/International Internship/Internship-25-26/scMINER_pipeline/RESULTS"

# Create a dedicated directory for results if it does not exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# Object to store validation scores across runs
silhouette_results <- list()

cat("--- Starting Automated MICA MDS Silhouette Analysis Pipelining ---\n")

# Loop through directories sequentially from k = 2 to k = 8
for (k in 2:8) {
  # Dynamically build the subfolder path on your OneDrive layout
  sub_folder <- file.path(base_mds_dir, paste0("mica_mds_output_all_genes_", k))
  
  # Automatically capture any generated clustering text file inside the folder
  txt_file <- list.files(sub_folder, pattern = "^clustering_UMAP.*\\.txt$", full.names = TRUE)[1]
  
  if (!is.na(txt_file) && file.exists(txt_file)) {
    # Read the tab-separated MICA text data
    mds_data <- read.table(txt_file, header = TRUE, sep = "\t")
    
    # Isolate the low-dimensional spatial coordinates and cluster integers
    coordinates_matrix <- as.matrix(mds_data[, c("X", "Y")])
    cluster_labels     <- as.integer(mds_data$label)
    
    # Calculate Euclidean distances and cell-level silhouette widths
    spatial_distances  <- dist(coordinates_matrix, method = "euclidean")
    sil_object         <- cluster::silhouette(x = cluster_labels, dist = spatial_distances)
    avg_score          <- mean(sil_object[, "sil_width"])
    
    # Store results for downstream comparison
    silhouette_results[[as.character(k)]] <- avg_score
    cat(paste0("Folder _", k, " -> Target File: ", basename(txt_file), " | Silhouette Score = ", round(avg_score, 4), "\n"))
  } else {
    cat(paste0("X Folder _", k, ": No valid clustering .txt file discovered.\n"))
  }
}

# Determine the mathematically optimal cluster configuration from the loop
if (length(silhouette_results) > 0) {
  best_k <- names(silhouette_results)[which.max(unlist(silhouette_results))]
  cat("\n==> Optimal configuration based on MDS parameters is K =", best_k, "with a validation score of", round(silhouette_results[[best_k]], 4), "!\n\n")
} else {
  stop("CRITICAL ERROR: No data files could be correctly accessed or scanned.")
}


# ==============================================================================
# STEP 6: DYNAMICALLY INTEGRATE WINNING MICA MDS LABELS & VISUALISE
# ==============================================================================
if (exists("best_k")) {
  winning_folder <- file.path(base_mds_dir, paste0("mica_mds_output_all_genes_", best_k))
  winning_file   <- list.files(winning_folder, pattern = "^clustering_UMAP.*\\.txt$", full.names = TRUE)[1]
  
  if (!is.na(winning_file) && file.exists(winning_file)) {
    winning_data <- read.table(winning_file, header = TRUE, sep = "\t")
    match_indices <- match(rownames(pData(combined.eset)), winning_data$ID)
    
    # Om jouw originele MICAplot code te hergebruiken, slaan we de clusters op als 'clusterID'
    pData(combined.eset)$clusterID <- as.factor(winning_data$label[match_indices])
    pData(combined.eset)$UMAP_1    <- winning_data$X[match_indices]
    pData(combined.eset)$UMAP_2    <- winning_data$Y[match_indices]
    pData(combined.eset)$sample    <- as.factor(pData(combined.eset)$sample)
    
    cat("--- SUCCESS: MICA OUTPUT DYNAMICALLY INTEGRATED (K =", best_k, ") ---\n")
    print(table(pData(combined.eset)$clusterID))
    
    # --------------------------------------------------------------------------
    # JOUW ORIGINELE VISUALISATIES (Nu volledig functioneel!)
    # --------------------------------------------------------------------------
    # 1. Plot gekleurd per Cluster ID
    print(MICAplot(input_eset = combined.eset, 
                   color_by = "clusterID", 
                   X = "UMAP_1", Y = "UMAP_2", 
                   point.size = 0.1, 
                   fontsize.cluster_label = 6))

    # 2. Plot gekleurd per Sample (Biofilm vs Planktonic) om batch effecten te checken
    print(MICAplot(input_eset = combined.eset, 
                   color_by = "sample", 
                   X = "UMAP_1", Y = "UMAP_2", 
                   point.size = 0.1))
    
  } else {
    stop(paste("ERROR: MICA output file not found in folder:", winning_folder))
  }
} else {
  stop("CRITICAL ERROR: 'best_k' is missing. Run the STEP 5 loop first!")
}


# ==============================================================================
# STEP 7: COMPUTE CLUSTER-SPECIFIC MARKER GENES & BUBBLEPLOT
# ==============================================================================
cat(paste0("Initiating differential expression profiling using scMINER on optimal clusters (K = ", best_k, ")...\n"))

# Wilcoxon-based differential expression analysis per cluster
de_res <- getDE(input_eset  = combined.eset, 
                group_by    = "clusterID", 
                use_method  = "wilcoxon") # Change to "limma" if preferred

# Select the top 10 marker genes per cluster, ranked by log2FC
cluster_markers <- getTopFeatures(input_table     = de_res, 
                                  number          = 10, 
                                  group_by        = "g1_tag", 
                                  sort_by         = "log2FC", 
                                  sort_decreasing = TRUE)

print(head(cluster_markers))

# Export results to the RESULTS directory using the dynamic K-value in the filename
output_file_path <- file.path(results_dir, paste0("MICA_MDS_Cluster_K", best_k, "_Top10_Markers.csv"))
write.csv(cluster_markers, file = output_file_path, row.names = FALSE)
cat(paste0("SUCCESS: Markers exported to: ", basename(output_file_path), "\n"))

# --------------------------------------------------------------------------
# BUBBLEPLOT VISUALISATION OF TOP MARKER GENES ACROSS CLUSTERS
# --------------------------------------------------------------------------

# 1. Generate the base bubble plot with vertical labels (90 degrees) to prevent overlaps
p <- feature_bubbleplot(input_eset   = combined.eset, 
                        features     = unique(cluster_markers$feature), 
                        group_by     = "clusterID", 
                        xlabel.angle = 90)

# 2. Fine-tune X-axis text for perfect PowerPoint legibility
p_fixed <- p + theme(
  axis.text.x = element_text(
    size = 8,          # Shrink text slightly so more genes fit side-by-side
    hjust = 1,         # Right-align labels so they line up perfectly under each dot column
    vjust = 0.5        # Center text vertically relative to the rotation point
  )
)

# 3. Display the corrected plot on screen
print(p_fixed)

# 4. Export high-res image using your relative path
# This automatically saves it into the 'PLOT' folder of your current working directory
ggsave(filename = "./PLOT/top_markers_widescreen.png", plot = p_fixed, width = 14, height = 6, dpi = 300)

# ==============================================================================
# STEP 8: CELLTYPIST INTEGRATION AND BIOLOGICAL CONSENSUS MAPPING
# ==============================================================================
if (file.exists(celltypist_file)) {
  celltypist_data <- read.table(celltypist_file, header = TRUE, sep = "\t")
  
  # Map high-confidence CellTypist annotations to the synchronized dataset cells
  pData(combined.eset)$Celltypist_Label <- as.factor(
    celltypist_data$predicted_labels[match(rownames(pData(combined.eset)), 
                                           celltypist_data$ID)]
  )
  
  # Compute the downstream consensus matrix (cross-tabulation evaluation)
  consensus_table <- table(pData(combined.eset)$clusterID,        # ← fixed
                           pData(combined.eset)$Celltypist_Label)
  
  cat("\n--- CONSENSUS ANNOTATION MATRIX (MICA MDS vs CellTypist) ---\n")
  print(consensus_table)
  cat("-------------------------------------------------------------\n")
  
  # Save the cross-tabulation table for the final report appendix
  write.csv(consensus_table, 
            file = file.path(results_dir, "MICA_vs_Celltypist_Consensus.csv"))
  
  # Generate publication-quality proportional distribution plot using ggplot2
  ggplot(pData(combined.eset), aes(x = clusterID,                 # ← fixed
                                   fill = Celltypist_Label)) +
    geom_bar(position = "fill") +
    theme_minimal() +
    labs(title = paste0("CellTypist Label Distribution Across MICA MDS Clusters (K=", best_k, ")"),
         x     = "MICA MDS Clusters",
         y     = "Proportion of Cells",
         fill  = "CellTypist Label") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  ggsave(file.path(results_dir, "MICA_Celltypist_Consensus_Plot.png"), 
         width = 8, height = 5)
  cat("Consensus visualization saved to the RESULTS folder!\n")
  
} else {
  cat("\n[Notice]: CellTypist prediction file was not found at:", celltypist_file, "\n")
  cat("Pipeline is ready; step 8 will execute once CellTypist output is available.\n")
}

markers <- read.csv("./RESULTS/MICA_MDS_Cluster_K4_Top10_Markers.csv")
View(markers)  # opent netjes in RStudio


# Figure 1 (CellTypist) defines the biological state of each mathematical cluster,
# while Figure 2 (Bubble Plot) uncovers the molecular mechanisms that drive those states.
#
# Starting with Cluster 1, the Step 8 bar chart confirms it is almost exclusively
# planktonic. The bubble plot then reveals why: both cntA and the nickel ABC transporter
# protein are uniquely and heavily expressed in this cluster, making cntA a robust,
# high-confidence biomarker for planktonic cells.
#
# Cluster 4 tells a different story. The bar chart shows it is heavily enriched with
# biofilm cells, and the bubble plot exposes the underlying machinery: rpoA (an RNA
# polymerase subunit) and 6S RNA (a regulatory RNA) are both highly active here.
# This points to a major transcriptional overhaul at the moment of biofilm commitment,
# likely orchestrated by 6S RNA-mediated regulation during the stress and switching phase.
#
# Taken together, the CSV data and bar chart validate the cell populations, while the
# bubble plot serves as the visual discovery tool that exposes which bacterial machinery
# is actively at work — a central finding of this international internship project.




