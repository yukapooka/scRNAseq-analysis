## Seurat Guided Clustering Tutorial

#####################################
## Setup the Seurat Object
#####################################


## Install and load packages
# Install.packages('Seurat')
#reticulate::py_install(packages = 'umap-learn')
library(dplyr)
library(Seurat)
library(patchwork)

## Load data
# Load PBMC data
pbmc.data<- Read10X(data.dir="data/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts=pbmc.data, 
                           project = "pbmc3k", 
                           min.cells = 3, 
                           min.features = 200)

# Visualize the genes in first 30 cells
pbmc.data[1:5,1:30]
dense.size <- object.size(as.matrix(pbmc.data))
dense.size
sparse.size <- object.size(pbmc.data)
sparse.size
dense.size/sparse.size


#####################################
## Pre-processing workflow
#####################################

## QC
# Percent reads mapped to mt genome
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# Show QC metrics for the first 5 cells
head(pbmc@meta.data,5)

# Visualize QC metrics as violin plot
VlnPlot(pbmc, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"),
         ncol = 3)
# FeatureScatter to visualize feature-feature relatonships
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter cells that have unique feature counts > 2,500 or < 200 or >5% mt counts
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

## Normalization by "LogNormalize" method, i.e. normalizes the feature expression 
## measurements for each cell by the total expression, multiplies this by 
## a scale factor (10,000 by default), and log-transforms

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

## Identify highly variable features (feature selection)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the top 10 variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# Plot with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot=plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scaling the data
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")

# Perform linear dimensional reduction (PCA) on the scaled data
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results in a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca")

# Plot for PCA 1
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

# Plot for PCAs 1-15
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

## Determine the 'dimensionality' of the dataset
# The top principal components therefore represent a robust compression of the dataset. 
# JackStraw approach: We randomly permute a subset of the data (1% by default) and 
# rerun PCA, constructing a ‘null distribution’ of feature scores, and repeat this procedure. 
# We identify ‘significant’ PCs as those who have a strong enrichment of low p-value features.

# Might take a while to run this.
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
#JackStrawPlot() function provides a visualization tool for comparing the distribution of p-values 
# for each PC with a uniform distribution (dashed line).
# ‘Significant’ PCs will show a strong enrichment of features with 
# low p-values (solid curve above the dashed line). In this case it appears that there is a 
# sharp drop-off in significance after the first 10-12 PCs.

# An alternative heuristic method generates an ‘Elbow plot’: a ranking of principle components based on 
# the percentage of variance explained by each one
ElbowPlot(pbmc)
# Majority of true signal is captured in the first 10 PC
# The first is more supervised, exploring PCs to determine relevant sources of heterogeneity, and could be 
# used in conjunction with GSEA for example. The second implements a statistical test based on a random null model, 
# but is time-consuming for large datasets, and may not return a clear PC cutoff. The third is a heuristic 
# that is commonly used, and can be calculated instantly. In this example, all three approaches yielded similar results, 
# but we might have been justified in choosing anything between PC 7-12 as a cutoff.

#####################################
## Cluster the cells
#####################################

# Graph-based clustering approach: embed cells in a graph structure - 
# for example a K-nearest neighbor (KNN) graph, with edges drawn between cells with similar feature expression patterns, 
# and then attempt to partition this graph into highly interconnected ‘quasi-cliques’ or ‘communities’.
# first construct a KNN graph based on the euclidean distance in PCA space, and refine the edge weights between any two cells 
# based on the shared overlap in their local neighborhoods (Jaccard similarity).

# To cluster the cells, we next apply modularity optimization techniques such as the Louvain algorithm (default) or SLM, to iteratively group cells together, 
# with the goal of optimizing the standard modularity function.

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

#Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

## Run non-linear Dim Reduction (UMAP, tSNE)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")
saveRDS(pbmc, file = "output/pbmc_tutorial.rds")

#####################################
## Finding DE features (cluster biomarkers)
#####################################
# Find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)

cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
# you can plot raw counts as well
VlnPlot(pbmc, features = c("MS4A1", "CD79A"), slot = "counts", log = TRUE)

FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
                               "CD8A"))

top10 <- pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
# DoHeatmap(pbmc, features = top10$gene) + NoLegend()
DoHeatmap(pbmc, features = top10$gene)

#####################################
## Assign cell type identity to clusters
#####################################
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
saveRDS(pbmc, file = "output/pbmc3k_final.rds")
Í›