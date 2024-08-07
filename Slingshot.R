library(R.utils)
library(Seurat)
library(tradeSeq)
library(dplyr)
library(slingshot)
library(Matrix)
library(patchwork)
library(tidyverse)


#Create data folder
dir.create("data")

# URL of the file to be downloaded
file_url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72857/suppl/GSE72857%5Fumitab%2Etxt%2Egz"

# Destination file path
dest_file <- "data/GSE72857_umitab.txt.gz"

# Download the file
download.file(url = file_url, destfile = dest_file, mode = "wb")

# Print a message to indicate the download is complete
cat("File downloaded to", dest_file, "\n")

# Decompress the file
gunzip("data/GSE72857_umitab.txt.gz", remove = FALSE)

# Loading matrix into R and converting it to compressed .rds format.

data <- read.delim("data/GSE72857_umitab.txt", header = T, row.names = 1)
comp_matrix <- Matrix::Matrix(as.matrix(data), sparse = T)
saveRDS(comp_matrix, "data/GSE72857_umitab.rds")


View(umi_counts)

# Loading data
umi_counts <- readRDS("data/GSE72857_umitab.rds")
umi_counts <- umi_counts[, c(T, F, F, F, F)]
dim(umi_counts)


# Load the raw accession meta data to merge required columns to our seurat object
cell_data <- read.csv("C:\\Users\\divya_vq9ublx\\Downloads\\SraRunTable.txt")
View(cell_data)

# Create the new column names and merge the data
data$Sample_name <- cell_data$Sample.Name
data$Cell_labels <- cell_data$source_name

View(data@meta.data)


# Define a color pallete to use
pal <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))

#.......................Seurat Workflow........................................

# Data analysis with Seurat pipeline
data <- CreateSeuratObject(counts = umi_counts)
data <- NormalizeData(data)
data <- FindVariableFeatures(data, nfeatures = 2000)
data <- ScaleData(data)
data <- RunPCA(data)
data <- FindNeighbors(data)
data <- FindClusters(data, resolution = 1)
data <- RunUMAP(data, n.neighbors = 10, dims = 1:50, spread = 2, min.dist = 0.3)



View(data@meta.data)

# Plot the clusters
DimPlot(data, reduction = 'umap', group.by = "Cell_labels")
DimPlot(data, reduction = 'umap', group.by = "seurat_clusters")

#Finding differentially expressed features (cluster biomarkers)
# find all markers of cluster 2
cluster2_markers <- FindMarkers(data, ident.1 = 2)
head(cluster2_markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5_markers <- FindMarkers(data, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5_markers, n = 5)


# find markers for every cluster compared to all remaining cells
data_markers <- FindAllMarkers(data, only.pos = TRUE)
data_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
cluster0_markers <- FindMarkers(data, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

data_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10_data
DoHeatmap(data, features = top10_data$gene) + NoLegend()

View(top10_data)

saveRDS(top10_data,"~/DataAnalysis/data/top10_data.rds")

View(data)


View(data@reductions$umap@cell.embeddings)

# Save the objects as separate matrices for input in slingshot
dimred <- data@reductions$umap@cell.embeddings
clustering <- Idents(data)
clustering <- data$RNA_snn_res.1
counts <- as.matrix(data@assays$RNA$counts[data@assays$RNA@meta.data$var.features])

View(counts)
# Get the variable features
var_features <- VariableFeatures(data)

# Extract the counts matrix for the variable features
counts <- as.matrix(data@assays$RNA$counts[var_features, ])

View(dimred)


#.............Trajectory inference with Slingshot...........................

# Defining cell lineages with Slingshot

suppressPackageStartupMessages({
  library(slingshot)
})
library(slingshot)

# Run default Slingshot lineage identification

lineages <- getLineages(data = dimred, clusterLabels = clustering)

lineages

# Plot the lineages
par(mfrow = c(1, 2))
plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)
for (i in levels(clustering)) {
  text(mean(dimred[clustering == i, 1]), mean(dimred[clustering == i, 2]), labels = i, font = 2)
}

# Second plot: scatter plot with lineages
plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)
lines(lineages, lwd = 3, col = "black")

lineages <- getLineages(data = dimred,
                        clusterLabels = clustering,
                        #end.clus = c("11","7","10","9","5"), #define how many branches/lineages to consider
                        start.clus = "0") #define where to start the trajectories

lineages

#Plot the lineages
par(mfrow=c(1,2))
plot(dimred[,1:2], col = pal[clustering],  cex=.5,pch = 16)
for(i in levels(clustering)){
  text( mean(dimred[clustering==i,1]),
        mean(dimred[clustering==i,2]), labels = i,font = 2) }

# Extract UMAP embeddings
dimred <- Embeddings(data, "umap")

# Extract clustering
clustering <- Idents(data)

# Set up colors for clusters
pal <- rainbow(length(unique(clustering)))
names(pal) <- levels(clustering)

# Run Slingshot
lineages <- slingshot(data = dimred, clusterLabels = clustering, reducedDim = 'UMAP')

# Set up the plotting area
par(mfrow = c(1, 2))

# First plot: scatter plot with clustering
plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 16, main = "Clusters")
for (i in unique(clustering)) {
  text(mean(dimred[clustering == i, 1]), mean(dimred[clustering == i, 2]), labels = i, font = 2)
}

# Second plot: scatter plot with lineages
plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 16, main = "Lineages")



# Extract and plot lineage curves
lineage_curves <- slingCurves(lineages)
for (curve in lineage_curves) {
  lines(curve$s[, 1], curve$s[, 2], lwd = 3, col = "black")
}

#....Finding differentially expressed genes......

BiocParallel::register(BiocParallel::SerialParam())
BiocManager::install("tradeSeq")

library(tradeSeq)

# Removing some genes to speed up the computations
filt_counts <- counts[rowSums(counts > 5) > ncol(counts)/100, ]
dim(filt_counts)

curves <- getCurves(lineages, approx_points = 300, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)
curves

sce <- fitGAM(counts = as.matrix(filt_counts), sds = curves)

plotGeneCount(curves, filt_counts, clusters = clustering, models = sce)

# Define function to plot
library(dplyr)
plot_differential_expression <- function(feature_id) {
  feature_id <- pseudotime_association %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id)
  cowplot::plot_grid(plotGeneCount(curves, filt_counts, gene = feature_id[1], clusters = clustering, models = sce) + ggplot2::theme(legend.position = "none"),
                     plotSmoothers(sce, as.matrix(counts), gene = feature_id[1]))
}



# Genes that change with pseudotime

pseudotime_association <- associationTest(sce)
pseudotime_association$fdr <- p.adjust(pseudotime_association$pvalue, method = "fdr")
pseudotime_association <- pseudotime_association[order(pseudotime_association$pvalue), ]
pseudotime_association$feature_id <- rownames(pseudotime_association)

feature_id <- pseudotime_association %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id)
plot_differential_expression(feature_id)

#Genes that change between two pseudotime points

pseudotime_start_end_association <- startVsEndTest(sce, pseudotimeValues = c(0, 1))
pseudotime_start_end_association$feature_id <- rownames(pseudotime_start_end_association)

feature_id <- pseudotime_start_end_association %>% filter(pvalue < 0.05) %>% top_n(1, waldStat) %>% pull(feature_id)

plot_differential_expression(feature_id)

#Genes that are different between lineages

different_end_association <- diffEndTest(sce)
different_end_association$feature_id <- rownames(different_end_association)
feature_id <- different_end_association %>% filter(pvalue < 0.05) %>% arrange(desc(waldStat)) %>% dplyr::slice(1) %>% pull(feature_id)
plot_differential_expression(feature_id)

branch_point_association <- earlyDETest(sce)
branch_point_association$feature_id <- rownames(branch_point_association)

feature_id <- branch_point_association %>% filter(pvalue < 0.05) %>% arrange(desc(waldStat)) %>% dplyr::slice(1) %>% pull(feature_id)
plot_differential_expression(feature_id)




