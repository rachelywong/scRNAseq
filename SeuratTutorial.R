# Seurat Tutorial

library(dplyr)
library(Seurat)

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "../data/pbmc3k/filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

## An object of class Seurat 
## 13714 features across 2700 samples within 1 assay 
## Active assay: RNA (13714 features)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

# Normalize the data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# the top one normalizes data to a set standard
# or you can normalize with this:
pbmc <- NormalizeData(pbmc)

# Identification of highly variable features
# high cell-to-cell variation in the dataset
# FindVariableFeatures function --> 2000 features per dataset
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
# ^^^ creates the identification scatter plot

# Scaling the data
# results are stored in --> pbmc[["RNA"]]@scale.data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

# Visualize dimensional reduction genes
# visualizes top genes associated with reduction components
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

# DimPlot
# graphs the output of a dimensional reduction technique (PCA by default)
# cells are colored by their identity class
DimPlot(pbmc, reduction = "pca")

# DimHeatmap
# allows for easy exploration of the primary sources of heterogeneity
# cells and features are ordered according to PCA scores
# draws a heatmap focusing on a principal component
# both cells and genes are sorted by their principal component scores
# allows for a nice visualization of sources of heterogeneity in the dataset
# heterogeneity --> quality or state of being diverse in character or content
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
# dims = 1:15
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# Dimensionality of dataset
# Jackstraw procedure 
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
# Visualization for comparing the distribution of p-values for each PC
# with a uniform distribution (dotted line)
JackStrawPlot(pbmc, dims = 1:15)

# An alternative heuristic method generates an "Elbow Plot"
# which is a ranking of principle components based on the percentage
# variance explained by each (ElbowPlot function)
# we can find where the majority of the true signal is captured (in which PCs)
ElbowPlot (pbmc)

# Cluster the cells
# FindNeighbours function --> takes the previously defined dimensionality
# of the dataset (ex// first 10 PCs) as input and constructs a graph
# based on the euclidean distance in PCA space and refines the edge weights
# between any two cells based on the shared overlap in their local 
# neighbourhoods (Jaccard similarity)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
# FindClusters function --> clusters the cells and applies optimization techniques
# (Louvain algorithm is default or SLM) 
# single cell datasets of around 3K cells works well with 0.4-1.2 resolution
# parameter that sets the "granularity" of the clustering
# optimal resolution often increases for larger datasets
pbmc <- FindClusters(pbmc, resolution = 0.5)
# Idents function --> finds clusters
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

# UMAP/tSNE to visualize dataset of clusters
# Install UMAP --> py.install(packages = 'umap-learn')
py.install(packages = 'umap-learn') # doesn't work
install.packages("umap")
library("umap", lib.loc="/Library/Frameworks/R.framework/Versions/3.6/Resources/library")
install.packages("umap-learn")

# Run UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)

# Create plot
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "umap")
# creates the UMAP plot

install.packages("PythonInR") #install python

# Restart R after installing packages
install.packages('reticulate')
library('reticulate')
reticulate::py_install(packages = 'umap-learn')

# Find differentially expressed features
# (cluster biomarkers)
# find all markers of cluster 1 using:
# find all markers of cluster 1
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
# download and load magrittr and dplyr packages first to use %>%
install.packages("magrittr")
install.packages("dplyr")
library (magrittr)
library (dplyr)
# make sure to restart R every time you install packages

# Test for differential expression 
cluster1.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

# Visualize marker expressions using VlnPlot
VlnPlot(pbmc, features = c("MS4A1", "CD79A")) #example given in tutorial
VlnPlot(pbmc, features = c("Neat1", "Dsg1a")) #features from cluster 5 (distinguish from 0 and 3)

# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
VlnPlot(pbmc, features = c("Neat1", "Dsg1a"), slot = "counts", log = TRUE)

# FeaturePlot
# featureplot can be used to visualize feature expressions on a tSNE or PCA plot
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
                               "CD8A")) #tutorial example
FeaturePlot(pbmc, features = c("Neat1", "Dsg1a")) # example to use

# DoHeatmap
# generates an expression heatmap for given cells and features
# plot the top 20 markers for each cluster using:
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

# Assign cell type identity to clusters
# can use canonical markers to easily match the unbiased clustering to known
# cell types
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
# may need to change cluster IDs, cluster IDs used do not match dataset
# cluster IDs from Seurat tutorial, not sure which IDs to use
 
