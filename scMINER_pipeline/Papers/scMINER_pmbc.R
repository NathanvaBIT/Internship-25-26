#C:/Users/natha/OneDrive/Bio-informatica 25-26'/International Internship

# Install and load the scMINER package from GitHub
devtools::install_github("jyyulab/scMINER")
library(scMINER)

# Load the built-in PBMC14k raw count dataset (17986 genes x 14000 cells)
data("pbmc14k_rawCount") #read in the raw counts files 


# Create a SparseEset object from the raw count matrix, including cell and gene metadata
pbmc14k_raw.eset <- createSparseEset(input_matrix = pbmc14k_rawCount, projectID = "PBMC14k", addMetaData = TRUE)
head(pData(pbmc14k_raw.eset)) #phenotype (cel metadata)
head(fData(pbmc14k_raw.eset)) #gene metadata

# Read the true cell type labels for PBMC14k, used later for validation of clustering results
true_label <- read.table(system.file("extdata/demo_pbmc14k/PBMC14k_trueLabel.txt.gz", package = "scMINER"), header = T, row.names = 1, sep = "\t", quote = "", stringsAsFactors = FALSE)
head(true_label)

#5.1 QC report 
#Generate an interactive HTML QC report to assess data quality before filtering
drawSparseEsetQC(input_eset = pbmc14k_raw.eset, output_html_file = "C:/Users/natha/OneDrive/Bio-informatica 25-26'/International Internship/PLOT/pbmc14k_rawCount.html", overwrite = TRUE)

#5.2.1 data filtration with auto mode
# Remove low-quality cells and genes automatically based on QC metrics (both genes and cells are filtered)
pbmc14k_filtered.eset <- filterSparseEset(pbmc14k_raw.eset, filter_mode = "auto", filter_type = "both")


#6.0 Data normalisation
# Normalize raw counts to log2-CPM (counts per million) to correct for differences in sequencing depth between cells
pbmc14k_log2cpm.eset <- normalizeSparseEset(pbmc14k_filtered.eset, scale_factor = 1000000, log_base = 2, log_pseudoCount = 1)

#6.1 Save SparseEset object
# Save the normalized SparseEset object as an RDS file so it can be reloaded without rerunning all previous steps
dir.create("C:/Users/natha/OneDrive/Bio-informatica 25-26'/International Internship/DATA", 
           recursive = TRUE)
saveRDS(pbmc14k_log2cpm.eset, file = "C:/Users/natha/OneDrive/Bio-informatica 25-26'/International Internship/DATA/pbmc14k_log2CPM.rds")

#7.1 Run MICA 
# Unsupervised methods are used as cross-validation when reference studies are unavailable or unreliable,
# as they do not depend on prior knowledge of markers or cell type labels defined by other authors.
# In this script we use MICA (Mutual Information-based Clustering Analysis) as our unsupervised method,
# which clusters cells based on their gene expression profiles without needing predefined cell type labels.
## generate MICA input in h5ad format: anndata package is needed

dir.create("C:/Users/natha/OneDrive/Bio-informatica 25-26'/International Internship/MICA", 
           recursive = TRUE)
generateMICAinput(input_eset = pbmc14k_log2cpm.eset, output_file = "C:/Users/natha/OneDrive/Bio-informatica 25-26'/International Internship/MICA/micaInput.txt", overwrite = TRUE)


#7.2 Graph Embedding > large dataset of 13605 cells so we use ge mode
# Run MICA clustering in graph embedding (ge) mode via the terminal (NOT in R) because of the large dataset size
# mica ge -i "C:/Users/natha/OneDrive/Bio-informatica 25-26'/International Internship/MICA/micaInput.h5ad" -o "C:/Users/natha/OneDrive/Bio-informatica 25-26'/International Internship/MICA/micaOutput" -minr 0.1 -maxr 9.0 -ss 0.05 -nw 4

###7.4 Integrate MICA outputs into SparseEset object
# Load MICA clustering results: each row is a cell with its UMAP coordinates (X, Y) and assigned cluster label
micaOutput <- read.table(system.file("extdata/demo_pbmc14k/MICA/clustering_UMAP_euclidean_20_2.05.txt", package = "scMINER"), header = TRUE, sep = "\t", quote = "", stringsAsFactors = F)
head(micaOutput)
# micaOutput contains the MICA clustering results, one row per cell
# ID    = cell barcode (unique identifier of the cell)
# X     = UMAP coordinate on the X-axis
# Y     = UMAP coordinate on the Y-axis
# label = cluster number assigned to that cell by MICA


# Add MICA clustering results (UMAP_1, UMAP_2, clusterID) to the SparseEset object
pbmc14k_log2cpm.eset <- addMICAoutput(pbmc14k_log2cpm.eset, mica_output_file = system.file("extdata/demo_pbmc14k/MICA/clustering_UMAP_euclidean_20_2.05.txt", package = "scMINER"), visual_method = "umap")
head(pData(pbmc14k_log2cpm.eset))

#Clustering results (3 new columns: UMAP_1, UMAP_2, clusterID) will be added to SparseExpressionsSet object
pbmc14k_log2cpm.eset <- addMICAoutput(pbmc14k_log2cpm.eset, mica_output_file = system.file("extdata/demo_pbmc14k/MICA/clustering_UMAP_euclidean_20_2.05.txt", package = "scMINER"), visual_method = "umap")
head(pData(pbmc14k_log2cpm.eset))

#7.5 Visualisation the MICA output 
#7.5.1 Color-coded by cluster labels
# At this stage cluster numbers are not yet linked to biological cell type names (see Chapter 8 for annotation)

library(ggplot2)
MICAplot(input_eset = pbmc14k_log2cpm.eset, color_by = "clusterID", X = "UMAP_1", Y = "UMAP_2", point.size = 0.1, fontsize.cluster_label = 6)


#7.5.2 Color-coded by true label of cell types
# Colour cells by their known true cell type label to compare with unsupervised clustering results
MICAplot(input_eset = pbmc14k_log2cpm.eset, color_by = "trueLabel", X = "UMAP_1", Y = "UMAP_2", point.size = 0.1, fontsize.cluster_label = 4)

###ERROR: object color_by is here not found 
#7.5.3 Color-coded by nUMI, for QC purpose
MICAplot(input_eset = pbmc14k_log2cpm.eset, color_by = "trueLabel", X = "UMAP_1", Y = "UMAP_2", point.size = 0.1, fontsize.cluster_label = 4)

#7.5.4 Color-coded by nFeature, for QC purpose
MICAplot(input_eset = pbmc14k_log2cpm.eset, color_by = "nFeature", do.logTransform = TRUE, point.size = 0.1)

#7.5.4 Color-coded by  pctMito values, for QC purpose
# Colour cells by number of detected genes to check for low-quality cells with few expressed genes
MICAplot(input_eset = pbmc14k_log2cpm.eset, color_by = "pctMito", do.logTransform = FALSE, point.size = 0.1)


#8. Cell type annotation: here we will link the biological name to the cluster number
#8.1 Supervised cell type annotation
# Load the signature table containing marker genes and their weights per cell type
# Signature score = weighted mean expression of marker genes for a given cell type
signature_table <- read.table(system.file("extdata/demo_pbmc14k/PBMC14k_signatureTable.txt", package = "scMINER"), header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
head(signature_table)

# Visualise signature scores per cluster in a bubble plot to identify which cluster corresponds to which cell types
draw_bubbleplot(input_eset = pbmc14k_log2cpm.eset, signature_table = signature_table, group_by = "clusterID")

#8.1.2 Annotate using individual marker genes
# Select 2 well-known marker genes per cell type (7 cell types) plus CD3D and CD4 for further annotation
genes_of_interest <-c("CD14", "LYZ", "GZMB", "NKG7", "CD19", "MS4A1", "CD8A", "CD8B", "SELL", "CCR7", "IL2RA", "FOXP3", "IL7R", "S100A4", "CD3D", "CD4")

#8.1.2.1 feature visualization: violin plot
# Violin plot shows the distribution of marker gene expression per cluster
feature_vlnplot(input_eset = pbmc14k_log2cpm.eset, features = genes_of_interest, group_by = "clusterID", ncol = 4)

#8.1.2.2 feature visualization: box plot
# Box plot shows median and spread of marker gene expression per clusters
feature_boxplot(input_eset = pbmc14k_log2cpm.eset, features = genes_of_interest, group_by = "clusterID", ncol = 4)

#8.1.2.3 feature visualization: scatter plot
# UMAP scatter plot colours each cell by marker gene expression to see spatial distribution across clusters
feature_scatterplot(input_eset = pbmc14k_log2cpm.eset, features = genes_of_interest, ncol = 4, location_x = "UMAP_1", location_y =  "UMAP_2", point.size = 0.5, legend.key_height = 0.3, legend.key_width = 0.2, fontsize.legend_title = 8, fontsize.legend_text = 6, fontsize.axis_title = 8, legend.position = "none")

#8.1.2.4 feature visualization: bubble plot
# Bubble plot combines expression level and percentage of expressing cells per cluster for each marker gene
feature_bubbleplot(input_eset = pbmc14k_log2cpm.eset, features = genes_of_interest, group_by = "clusterID", xlabel.angle = 45)

#8.1.2.5 feature visualization: heatmap
## Heatmap of marker genes across clusters
feature_heatmap(input_eset = pbmc14k_log2cpm.eset, features = genes_of_interest, group_by = "clusterID", scale_method = "none", annotation_columns = c("trueLabel"))

##8.2 Unsupervised cell type annotation
# Unsupervised methods are used as cross-validation when reference studies are unavailable or unreliable,
# as they do not depend on prior knowledge of markers or cell type labels defined by other authors.
## 1. To perform differential expression analysis in a 1-vs-rest manner for all groups
de_res1 <- getDE(input_eset = pbmc14k_log2cpm.eset[500,], group_by = "clusterID", use_method = "limma")

head(de_res1)

## 2. To perform differential expression analysis in a 1-vs-rest manner for one specific group
de_res2 <- getDE(input_eset = pbmc14k_log2cpm.eset, group_by = "clusterID", g1 = c("1"), use_method = "limma")

## 3. To perform differential expression analysis in a rest-vs-1 manner for one specific group
de_res3 <- getDE(input_eset = pbmc14k_log2cpm.eset, group_by = "clusterID", g0 = c("1"), use_method = "limma")

## 4. To perform differential expression analysis in a 1-vs-1 manner for any two groups
de_res4 <- getDE(input_eset = pbmc14k_log2cpm.eset, group_by = "clusterID", g1 = c("1"), g0 = c("3"), use_method = "limma")

cluster_markers <- getTopFeatures(input_table = de_res1, number = 10, group_by = "g1_tag", sort_by = "log2FC", sort_decreasing = TRUE)
dim(cluster_markers)

head(cluster_markers)

celltype_map <- c(`1`="CD4TN", `2`="CD4TCM", `3`="CD8TN", `4`="NK", `5`="B", `6`="Monocyte", `7`="CD4Treg")
pbmc14k_log2cpm.eset$cell_type <- as.character(celltype_map[pbmc14k_log2cpm.eset$clusterID])
head(pData(pbmc14k_log2cpm.eset))

## Violin plot of marker genes across clusters
draw_barplot(input_eset = pbmc14k_log2cpm.eset, group_by = "cell_type", color_by = "trueLabel_full", xlabel.angle = 45)

saveRDS(pbmc14k_log2cpm.eset, file = "/your-path/PBMC14k/DATA/pbmc14k_log2CPM_annotated.rds")
