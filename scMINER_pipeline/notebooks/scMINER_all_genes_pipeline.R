"C:/Users/natha/OneDrive/Bio-informatica 25-26'/International Internship/scMINER"
# Install and load the scMINER package from GitHub + we downloaded already the scMINER from github.
#devtools::install_github("jyyulab/scMINER")
library(scMINER)
library(Matrix)
library(ggplot2)
library(anndata)

# Since both datasets (BF and P) were already filtered, quality-controlled and 
# normalised in Python (log1p-CPM values confirmed), and already concatenated  
# into a single combined dataset (adata_bacteria), we skip the QC report, 
# filtration, normalisation and concatenation steps that are part of the standard 
# scMINER pipeline. We directly load the combined dataset, ready for MICA clustering.



# STEP 1: Load the combined dataset
# readInput_h5ad() returns an AnnData object with:
#   $X   → expression matrix (cells x genes)
#   $obs → cell metadata (sample, leiden, cell_type)
#   $var → gene metadata

adata_combined <- readInput_h5ad(h5ad_file = "./DATA/bacteria_all_genes_combined.h5ad")

# Prefix toevoegen op basis van sample kolom
sample_labels <- adata_combined$obs$sample  # "bf" of "p"
cell_numbers  <- adata_combined$obs_names   # "68", "96", etc.

adata_combined$obs_names <- paste0(sample_labels, "_", cell_numbers)

# Verify
print(head(adata_combined$obs_names, 5))  # Moet bf_68, bf_96... tonen

# Verify dimensions after transpose — expected: 1015 genes x 7884 cells
print(dim(adata_combined))


# Extract and transpose expression matrix (cells x genes → genes x cells)
combined_mtx <- as(t(adata_combined$X), "dgCMatrix")
rownames(combined_mtx) <- rownames(adata_combined$var)  # ← this restores the bf_ / p_ gene names
colnames(combined_mtx) <- rownames(adata_combined$obs)  # ← this restores the cell namesµ

print(dim(combined_mtx))  # expected: 1015 genes x 7884 cel

# STEP 2: Create SparseEset object
# NOTE: the argument is 'cellData', NOT 'meta.data' (confirmed in documentation)
cell_metadata <- as.data.frame(adata_combined$obs)

combined.eset <- createSparseEset(input_matrix = combined_mtx,
                                  cellData     = cell_metadata,  # ← correct argument
                                  projectID    = "BaSSSh_bacteria",
                                  addMetaData  = TRUE)

dim(combined.eset)          # expected: 1015 x 7884
head(pData(combined.eset))  # verify sample, cell_type columns are present

# STEP 3: Export as MICA input file
generateMICAinput(input_eset = combined.eset,
                  output_file = "./MICA/MICA_all_genes_Input.txt")

  
# STEP 4: Run MICA in terminal (NOT in R — heavy computation!)
# To be executed in the terminal once MICA is installed:
# mica ge -i "../MICA/MICA_all_genes_Input.txt" -o "../MICA/micaOutput" -nw 6 -nc 6




# STEP 5: Read MICA output back into R
# Replace "x.xx" with the actual resolution from the MICA output filename
combined.eset <- addMICAoutput(combined.eset,
                               mica_output_file = "../MICA/micaOutput/clustering_UMAP_euclidean_20_x.xx.txt",
                               visual_method = "umap")
head(pData(combined.eset))

# STEP 6: Visualise clustering results
# Color by cluster ID
MICAplot(input_eset = combined.eset, 
         color_by = "clusterID", 
         X = "UMAP_1", Y = "UMAP_2", 
         point.size = 0.1, 
         fontsize.cluster_label = 6)

# Color by sample (BF vs P) to check batch mixing
MICAplot(input_eset = combined.eset, 
         color_by = "sample", 
         X = "UMAP_1", Y = "UMAP_2", 
         point.size = 0.1)

# STEP 7: Unsupervised cell type annotation via differential expression
# Since no ground truth labels are available for the bacteria dataset,
# we use differential expression analysis to identify cluster-specific marker genes.
de_res <- getDE(input_eset = combined.eset, 
                group_by = "clusterID", 
                use_method = "limma")


# Get top 10 marker genes per cluster
cluster_markers <- getTopFeatures(input_table = de_res, 
                                  number = 10, 
                                  group_by = "g1_tag", 
                                  sort_by = "log2FC", 
                                  sort_decreasing = TRUE)

head(cluster_markers)

# Visualise marker genes across clusters
feature_bubbleplot(input_eset = combined.eset, 
                   features = unique(cluster_markers$feature), 
                   group_by = "clusterID", 
                   xlabel.angle = 45)

# STEP 8: Network inference (SJARACNe) — optional, to be run after clustering validation
# This step is computationally very intensive and should be run on a powerful machine.
# generateSJARACNeInput(combined.eset, output_dir = "../SJARACNe/")