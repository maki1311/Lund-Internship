devtools::install_github("satijalab/seurat", ref = "spatial")

library(Seurat)
library(SeuratData)
library(ggplot2)
library(cowplot)
library(dplyr)
library(cellexalvrR)

#download spatial data
InstallData("stxBrain")
brain <- LoadData("stxBrain", type = "anterior1")

#look for heterogenity
plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
plot_grid(plot1, plot2)

#data normalization
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)

#visualization of gene expression
SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))
p1 <- SpatialFeaturePlot(brain, features = "Ttr", pt.size.factor = 1)
p2 <- SpatialFeaturePlot(brain, features = "Ttr", alpha = c(0.1, 1))
plot_grid(p1, p2)

#dimensional reduction, clustering and visualization
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)
p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
plot_grid(p1, p2)
SpatialDimPlot(brain, do.hover = TRUE)
SpatialDimPlot(brain, cells.highlight = CellsByIdentities(object = brain, idents = c(1, 2, 5, 3, 
                                                                                     4, 8)), facet.highlight = TRUE, ncol = 3)
LinkedDimPlot(brain)

#Identification of spacially variable features
de_markers <- FindMarkers(brain, ident.1 = 4, ident.2 = 6)
SpatialFeaturePlot(object = brain, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3)
brain <- FindSpatiallyVariableFeatures(brain, assay = "SCT", features = VariableFeatures(brain)[1:1000], 
                                       selection.method = "markvariogram")
top.features <- head(SpatiallyVariableFeatures(brain, selection.method = "markvariogram"), 6)
SpatialFeaturePlot(brain, features = top.features, ncol = 3, alpha = c(0.1, 1))

#Subset out anatomical regions
cortex <- subset(brain, idents = c(1, 2, 3, 5, 6, 7))
# now remove additional cells, use SpatialDimPlots to visualize what to remove
# SpatialDimPlot(cortex,cells.highlight = WhichCells(cortex, expression = image_imagerow > 400 |
# image_imagecol < 150))
cortex <- subset(cortex, image_imagerow > 400 | image_imagecol < 150, invert = TRUE)
cortex <- subset(cortex, image_imagerow > 275 & image_imagecol > 370, invert = TRUE)
cortex <- subset(cortex, image_imagerow > 250 & image_imagecol > 440, invert = TRUE)
p1 <- SpatialDimPlot(cortex, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(cortex, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
plot_grid(p1, p2)

#Integration with single cell data
allen_reference <- readRDS("/home/butlera/Projects/muir/seurat_objects/allen_brain.rds")
library(dplyr)
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% 
  RunUMAP(dims = 1:30)
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE) %>% RunPCA(verbose = FALSE)
DimPlot(allen_reference, group.by = "subclass", label = TRUE)
anchors <- FindTransferAnchors(reference = allen_reference, query = cortex, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$subclass, prediction.assay = TRUE, 
                                  weight.reduction = cortex[["pca"]])
cortex[["predictions"]] <- predictions.assay
DefaultAssay(cortex) <- "predictions"
SpatialFeaturePlot(cortex, features = c("L2/3 IT", "L4"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)
cortex <- FindSpatiallyVariableFeatures(cortex, assay = "predictions", features = rownames(cortex), 
                                        r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(cortex), 4)
SpatialPlot(object = cortex, features = top.clusters, ncol = 2)
SpatialFeaturePlot(cortex, features = c("Astro", "L2/3 IT", "L4", "L5 PT", "L5 IT", "L6 CT", "L6 IT", 
                                        "L6b", "Oligo"), pt.size.factor = 1, ncol = 2, crop = FALSE, alpha = c(0.1, 1))

#multiple slices 
brain2 <- LoadData("stxBrain", type = "posterior1")
brain2 <- SCTransform(brain2, assay = "Spatial", verbose = FALSE)
brain.merge1 <- merge(brain1, brain2)
DefaultAssay(brain.merge1) <- "SCT"
VariableFeatures(brain.merge1) <- c(VariableFeatures(brain1), VariableFeatures(brain2))
brain.merge1 <- RunPCA(brain.merge1, verbose = FALSE)
brain.merge1 <- FindNeighbors(brain.merge1, dims = 1:30)
brain.merge1 <- FindClusters(brain.merge1, verbose = FALSE)
brain.merge1 <- RunUMAP(brain.merge1, dims = 1:30)
DimPlot(brain.merge1, reduction = "umap", group.by = c("ident", "orig.ident"))
SpatialDimPlot(brain.merge1)
SpatialFeaturePlot(brain.merge1, features = c("Hpca", "Plp1"))


#load in all four slices and merge
brain1 <- LoadData("stxBrain", type = "anterior1")
brain1 <- SCTransform(brain1, assay = "Spatial", verbose = FALSE)

brain2 <- LoadData("stxBrain", type = "posterior1")
brain2 <- SCTransform(brain2, assay = "Spatial", verbose = FALSE)

brain3 <- LoadData("stxBrain", type = "anterior2")
brain3 <- SCTransform(brain3, assay = "Spatial", verbose = FALSE)

brain4 <- LoadData("stxBrain", type = "posterior2")
brain4 <- SCTransform(brain4, assay = "Spatial", verbose = FALSE)

brain.merge <- merge(brain1, y = c(brain2, brain3, brain4), add.cell.ids = c("A1", "P1", "A2", "P2"))
DefaultAssay(brain.merge) <- "SCT"

#change merged data
VariableFeatures(brain.merge) <- c(VariableFeatures(brain1), VariableFeatures(brain2), VariableFeatures(brain3), VariableFeatures(brain4))
brain.merge <- RunPCA(brain.merge, verbose = FALSE)
brain.merge <- FindNeighbors(brain.merge, dims = 1:30)
brain.merge <- FindClusters(brain.merge, verbose = FALSE)
brain.merge <- RunUMAP(brain.merge, dims = 1:30, n.components=3)


DimPlot(brain.merge, reduction = "umap", group.by = c("ident", "orig.ident"))
SpatialFeaturePlot(brain.merge, features = c("Penk"))

#make coordinates list
a1 <- cbind(brain.merge@images$anterior1@coordinates$col,brain.merge@images$anterior1@coordinates$row,1)
a2 <- cbind(brain.merge@images$anterior2@coordinates$col,brain.merge@images$anterior2@coordinates$row,2)
p1 <- cbind((brain.merge@images$posterior1@coordinates$col + 130),brain.merge@images$posterior1@coordinates$row,3)
p2 <- cbind((brain.merge@images$posterior2@coordinates$col + 130),brain.merge@images$posterior2@coordinates$row,4)
coords <- rbind(a1, a2, p1, p2)

#export as cellexalVR Object
brain.data <- GetAssayData(object = brain.merge) #Extract the expression data
drl <- list(UMAP=Embeddings(object = brain.merge, reduction = "umap"),slice = coords) # Put the UMAP coordinated into a list
meta <- make.cell.meta.from.df(brain.merge[[]], c("Cluster", "orig.ident")) # Make the metadata using just the "Tissue" column

cvr <- new("cellexalvrR",data=brain.data,drc=drl) # Initialise a new cellexalvrR object with the expression data and UMAP
cvr <- set.specie(cvr,"mouse") # Set the specie to Mouse
cvr <- addCellMeta2cellexalvr(cvr,meta) # Add the metadata to the cellexalvrR object
export2cellexalvr(cvr,"10x_spatial_2")


