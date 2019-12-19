library(cellexalvrR)
library(data.table)
library(Matrix)
library(Seurat)
library(umap)
library(rgl)
library(misc3d)

#load data
heart <- read.table("filtered_spots_count_matrix_all_weeks.tsv",sep="\t",header=T)
genes <- heart$X
heart <- heart[,-1]
rownames(heart) <- genes

stage <- c(rep(1,238), rep(2,1515), rep(3,1358))

#make Seurat
Seurat_heart <- CreateSeuratObject(heart)
Seurat_heart[["stage"]] <- stage
Seurat_heart <- NormalizeData(Seurat_heart, normalization.method = "LogNormalize")
Seurat_heart <- FindVariableFeatures(object = Seurat_heart, mean.function = ExpMean, dispersion.function = LogVMR, mean.cutoff = c(0.1, 10), dispersion.cutoff = c(0.5, Inf))
Seurat_heart <- ScaleData(Seurat_heart, vars.to.regress = "stage")
Seurat_heart <- RunICA(Seurat_heart, verbose = FALSE)
Seurat_heart <- RunPCA(Seurat_heart)
Seurat_heart <- FindNeighbors(Seurat_heart, reduction = "ica")
Seurat_heart <- FindClusters(object = Seurat_heart)
Seurat_heart <- RunTSNE(Seurat_heart, reduction = "ica", dims = 1:30, dim.embed = 3)
Seurat_heart <- RunUMAP(Seurat_heart, dims = 1:30, n.components=3)

#plot Seurat
DimPlot(object = Seurat_heart, reduction = "umap", group.by = "seurat_clusters", pt.size = 1)
heatmap(Seurat_heart@reductions$ica)

#make coordinates
coords.list = NULL

cells <- names(heart)

z <- c("numeric", length = nrow(cells))
x <- c("numeric", length = nrow(cells))
y <- c("numeric", length = nrow(cells))

for (i in 1:length(cells)){
  z[i] <- as.numeric(strsplit(cells, "x") [[i]]) [1]
  x[i] <- as.numeric(strsplit(cells, "x") [[i]]) [2]
  y[i] <- as.numeric(strsplit(cells, "x") [[i]]) [3]
  
  coords.list <- cbind(as.numeric(x), as.numeric(y), as.numeric(z))
  rownames(coords.list) = NULL
  
}

#export as cellexalVR Object
heart.data <- GetAssayData(object = Seurat_heart) #Extract the expression data
drl <- list(UMAP=Embeddings(object = Seurat_heart, reduction = "umap"), TSNE=Embeddings(object = Seurat_heart, reduction = "tsne"),slice = coords.list) # Put the UMAP coordinated into a list
meta <- make.cell.meta.from.df(Seurat_heart[[]], c("seurat_clusters")) # Make the metadata using just the "Tissue" column

cvr <- new("cellexalvrR",data=heart.data,drc=drl) # Initialise a new cellexalvrR object with the expression data and UMAP
cvr <- set.specie(cvr,"human") # Set the species
cvr <- addCellMeta2cellexalvr(cvr,meta) # Add the metadata to the cellexalvrR object
export2cellexalvr(cvr,"heart_spatial")
