### this is the pipeline test, we're using it on the pmbc14k_integrated dataset 
"C:/Users/natha/OneDrive/Bio-informatica 25-26'/International Internship/scMINER"
# Install and load the scMINER package from GitHub + we downloaded already the scMINER from github.
#devtools::install_github("jyyulab/scMINER")
library(scMINER)
library(Matrix)
library(ggplot2)
library(anndata)

# STEP 1: Load the built-in PBMC14k dataset
# This dataset contains 17986 genes x 14000 cells, with 7 known cell types:
# B, CD4TN, CD4TCM, CD4Treg, CD8TN, NK, Monocyte
data("pbmc14k_rawCount")
pbmc <- pbmc14k_rawCount
dim(pbmc)

## read the true labels of cell type for PBMC14k dataset
true_label <- read.table(system.file("extdata/demo_pbmc14k/PBMC14k_trueLabel.txt.gz", package = "scMINER"), 
                         header = TRUE, row.names = 1, sep = "\t", 
                         quote = "", stringsAsFactors = FALSE)
head(true_label)

# STEP 3: Create SparseEset object
pbmc_eset <- createSparseEset(input_matrix = pbmc,
                              cellData     = true_label,
                              projectID    = "PBMC14k",
                              addMetaData  = TRUE)

dim(pbmc_eset)          # expected: 17986 x 14000
head(pData(pbmc_eset))  # verify trueLabel column is present


# STEP 4: Data normalization
pbmc_eset <- normalizeSparseEset(pbmc_eset)
# STEP 5: Generate MICA input file
# MICA will cluster cells based on mutual information distances
generateMICAinput(input_eset  = pbmc_eset,
                  output_file = "./MICA/PBMC14k_MICA_Input.txt")





##---------------------------------------------------------------##
# STEP 6: Run MICA in terminal (NOT in R — heavy computation!)
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