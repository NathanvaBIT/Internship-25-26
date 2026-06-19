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
  best_k <- "7"
  cat("\n==> Manually set configuration: K =", best_k, "\n\n")
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

# ==============================================================================
# UMAP PLOTS PER MARKER GENE
# ==============================================================================

library(viridis)

# --- Project Directory Configuration ---
plot_dir <- "C:/Users/natha/OneDrive/Bio-informatica_25-26/International Internship/Internship-25-26/scMINER_pipeline/PLOT/"
dir.create(plot_dir, showWarnings = FALSE)

# Use combined.eset directly — DO NOT reload .RData (would overwrite clusterID!)
active_eset <- combined.eset

# Pull UMAP coordinates from pData (already integrated in Step 6)
umap_df <- data.frame(
  UMAP_1 = pData(active_eset)$UMAP_1,
  UMAP_2 = pData(active_eset)$UMAP_2
)

# Extract the expression matrix
expr_matrix <- exprs(active_eset)

# Define the list of 53 custom marker genes
marker_genes <- c(
  "nickel ABC transporter substrate-binding protein", "accD", "nusA", "6S RNA", 
  "SAUSA300_RS07045", "SAUSA300_RS08485", "pxpA", "tpiA", "rpoY", "pyrH", 
  "thiW", "tcyP", "sasA", "mvaK2", "sodium:proton antiporter-1", "clpB", 
  "HAD family hydrolase", "MFS transporter-3", "dapA", "SAUSA300_RS02370", 
  "alanine glycine permease-1", "acyltransferase", 
  "bifunctional metallophosphatase/5'-nucleotidase-1", "SAUSA300_RS09795", 
  "int", "sek", "icaB", "queE", "lpl9", "tandem-type lipoprotein-1", 
  "ATPase-1", "N-acetyltransferase-2", "haloacid dehalogenase-2", 
  "energy coupling factor transporter S component ThiW", "rpoA", "glpT", 
  "gpmA", "pbp1", "polX", "cls2", "accA", "ebh", "cntA", 
  "signal recognition particle sRNA large type", "tuf", "nasE", 
  "SAUSA300_RS15635", "xylose isomerase", "tagO", "tmaH", "pcrB", "pepQ2"
)

# Cross-reference with the active dataset
valid_genes <- marker_genes[marker_genes %in% rownames(active_eset)]
print(paste("Validated:", length(valid_genes), "out of 53 genes found in active dataset."))

# Loop through verified marker genes to create separate PDFs
print("Generating individual viridis-highlighted PDFs...")

for (gene in valid_genes) {
  
  # Create Windows-safe file names
  safe_gene_name <- gsub("[^A-Za-z0-9_-]", "_", gene)
  gene_pdf_path  <- paste0(plot_dir, "MICA_TOP_MARKER_GENES_Leiden_clustering_", safe_gene_name, ".pdf")
  
  umap_df$Expression <- expr_matrix[gene, ]
  
  pdf(gene_pdf_path, width = 8, height = 7)
  
  p <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = Expression)) +
    geom_point(size = 0.6, alpha = 0.8) +
    scale_color_viridis(option = "D", direction = 1) + 
    labs(title = paste("Expression of", gene),
         x = "UMAP_1", y = "UMAP_2",
         color = "Log Expression") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      panel.grid = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank()
    )
  
  print(p)
  dev.off()
}

print(paste("\n--- SUCCESS: All", length(valid_genes), "individual PDFs generated in:", plot_dir, "---"))

# ==============================================================================
# BUBBLEPLOT HIGHLY EXPRESSED COMMON GENES MICA
# ==============================================================================

# Load the highly expressed common genes from txt file
common_genes_path <- "C:/Users/natha/OneDrive/Bio-informatica_25-26/International Internship/Internship-25-26/scMINER_pipeline/RESULTS/highly_expressed_genes_common.txt"

common_genes <- readLines(common_genes_path)
common_genes <- trimws(common_genes[common_genes != ""])
cat(paste0("Loaded ", length(common_genes), " genes from file.\n"))

# Cross-reference with active dataset
valid_common_genes <- common_genes[common_genes %in% rownames(active_eset)]
cat(paste0("Validated: ", length(valid_common_genes), " out of ", 
           length(common_genes), " genes found in active dataset.\n"))

if (length(valid_common_genes) == 0) {
  stop("ERROR: No genes from txt file found in dataset. Check gene names!")
}

# Compute average expression and percentage of expressing cells per cluster
bubble_data <- do.call(rbind, lapply(levels(pData(active_eset)$clusterID), function(cl) {
  idx     <- pData(active_eset)$clusterID == cl
  sub_mat <- expr_matrix[valid_common_genes, idx, drop = FALSE]
  data.frame(
    gene    = valid_common_genes,
    cluster = cl,
    avg_exp = rowMeans(sub_mat),
    pct_exp = rowMeans(sub_mat > 0) * 100
  )
}))

# Order clusters numerically on x-axis
bubble_data$cluster <- factor(bubble_data$cluster,
                              levels = sort(unique(as.numeric(bubble_data$cluster))))

# Generate the bubbleplot
p_bubble <- ggplot(bubble_data, aes(x = cluster, y = gene, size = pct_exp, color = avg_exp)) +
  geom_point() +
  scale_color_viridis_c(option = "D", name = "Avg Expression") +
  scale_size(range = c(1, 10), name = "% Expressing Cells") +
  labs(
    title    = "Bubbleplot Duplicated Genes MICA",
    subtitle = paste0("Highly expressed common genes across ",
                      nlevels(pData(active_eset)$clusterID),
                      " clusters (n = ", length(valid_common_genes), " genes)"),
    x        = "Cluster",
    y        = "Gene"
  ) +
  theme_minimal() +
  theme(
    plot.title      = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle   = element_text(size = 10, hjust = 0.5, color = "grey40"),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y     = element_text(size = 9),
    panel.grid      = element_line(color = "grey90"),
    legend.position = "right"
  )

# Export bubbleplot to PDF in RESULTS folder
bubble_pdf_path <- file.path(results_dir, "Bubbleplot_Duplicated_Genes_MICA.pdf")
pdf(bubble_pdf_path, width = 10, height = 7)
print(p_bubble)
dev.off()
cat(paste0("SUCCESS: Bubbleplot exported to: ", basename(bubble_pdf_path), "\n"))

# ==============================================================================
# SAVE FINAL MICA ESET OBJECT (with clusterID)
# ==============================================================================
saveRDS(combined.eset,
        file = file.path(results_dir, "MICA_combined_eset_K7.rds"))
cat("SUCCESS: MICA eset saved as MICA_combined_eset_K7.rds\n")

# ==============================================================================
# ==============================================================================
# ==============================================================================

# 1. Zet je werkmap correct naar jouw projectmap


# 2. Laad de benodigde bibliotheek en je opgeslagen object uit de RESULTS map
library(Biobase)
combined.eset <- readRDS("./RESULTS/MICA_combined_eset_K7.rds")

# 3. Trek de getransponeerde expressiematrix eruit (Cellen x Genen)
matrix_transposed <- t(exprs(combined.eset))

# 4. Trek de cell metadata (pData) eruit
metadata_table <- pData(combined.eset)

# Waarom: Veiligheidscheck om er zeker van te zijn dat de coördinaten erin zitten
if (!"UMAP_1" %in% colnames(metadata_table)) {
  stop("CRITICAL ERROR: UMAP_1 (MICA coördinaten) niet gevonden in pData! Run eerst Step 6 van je script.")
}

# 5. Schrijf de data weg als CSV naar je huidige map
write.csv(as.matrix(matrix_transposed), file="./matrix_voor_python.csv", row.names=TRUE)
write.csv(metadata_table, file="./metadata_voor_python.csv", row.names=TRUE)

print("SUCCESS: Deel 1 afgerond! De CSV-bestanden staan klaar voor Python.")



import os
import anndata as ad
import pandas as pd
import numpy as np

# 1. Zorg dat Python in exact dezelfde map werkt
current_dir = "C:/Users/natha/OneDrive/Bio-informatica_25-26/International Internship/Internship-25-26/scMINER_pipeline"
os.chdir(current_dir)

# 2. Lees de zojuist gemaakte CSV's in
print("Data inladen...")
counts = pd.read_csv("./matrix_voor_python.csv", index_col=0)
metadata = pd.read_csv("./metadata_voor_python.csv", index_col=0)

# Waarom: Sanity check. De rijen van de matrix moeten 100% matchen met de rijen van de metadata.
if not (counts.index == metadata.index).all():
  raise ValueError("FOUT: De Cel-ID's van de matrix en metadata komen niet overeen!")

# 3. Het AnnData object initialiseren
# Waarom: We stoppen de counts in 'X' en koppelen de metadata aan 'obs'
adata = ad.AnnData(X=counts.values, obs=metadata)

# Waarom: Expliciet de celnamen en gennamen toewijzen zodat indexering werkt
adata.obs_names = counts.index.astype(str)
adata.var_names = counts.columns.astype(str)

# 4. HIER GEBEURT DE FIX: MICA-coördinaten toewijzen aan de UMAP embedding matrix
# Waarom: We trekken de UMAP_1 en UMAP_2 kolommen los uit adata.obs, 
# zetten ze om naar een numpy matrix, en slaan ze op in .obsm['X_umap']
adata.obsm['X_umap'] = adata.obs[['UMAP_1', 'UMAP_2']].to_numpy()

# Waarom: Nu de coördinaten veilig in .obsm['X_umap'] staan, verwijderen we de 
# losse kolommen uit .obs om je metadata tabel clean en overzichtelijk te houden.
adata.obs = adata.obs.drop(columns=['UMAP_1', 'UMAP_2'])

# 5. Opslaan als .h5ad bestand
# Waarom: Dit is de finale stap. We slaan het object op in de DATA map.
os.makedirs("./DATA", exist_ok=True)
output_file = "./DATA/bacteria_all_genes_with_mica_coordinates.h5ad"

adata.write_h5ad(output_file)

print(f"\n--- SUCCESS ---")
print(f"Je AnnData object is succesvol geconverteerd en opgeslagen!")
print(f"Bestandslocatie: {output_file}")
print(f"Data structuur: {adata.shape[0]} cellen x {adata.shape[1]} genen")
print(f"Beschikbare coördinaten: {list(adata.obsm.keys())} (bevat je MICA data)")
print(f"Beschikbare metadata: {adata.obs.columns.tolist()} (bevat o.a. 'clusterID' en 'sample')")

