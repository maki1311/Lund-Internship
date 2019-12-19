library("rgl")
library("misc3d")
library("ggplot2")
library("dplyr")

#gridpoints for shaping the brain and make text file for VR
cubegrid <- expand.grid(1:50,1:50,1:50)
cubegrid_data <- read.table("cubegrid.txt", header = TRUE, row.names=1)

#import gridpoints for the brain selected in VR
raw_brain_grid <- read.delim("selection_rawbrain.txt", header = FALSE, sep = "\t")

#select the datapoints inside the brain
index_raw_brain <- raw_brain_grid$V1[which(raw_brain_grid$V4 == 0)]

#use index for first_brain to select the correct lines in original raw brain grid
#make column containing the CellID Information
ID <- c(1:125000)
CellID_raw <- ID[which(raw_brain_grid$V4 == 0)]

#make matrix that contains the brain points and the CellID for VR
first_brain <- matrix(data = c(CellID_raw,cubegrid_data$Dim1[index_raw_brain],cubegrid_data$Dim2[index_raw_brain],cubegrid_data$Dim3[index_raw_brain])
                      , nrow = 28600, ncol = 4, byrow = FALSE)
colnames(first_brain) <- c("CellID", "x", "y", "z")

#make txt file containing the new matrix and reshape the brain in VR
write.table(first_brain, "raw_braingrid", sep = "\t")

#import reshaped brain txt file to VR and make further selections resulting in the final brain grid
#reshape brain
reshaped_brain <- read.delim("selection_reshapedbrain.txt", header = FALSE, sep = "\t")
brain_index <- reshaped_brain$V1[which(reshaped_brain$V4 == 1)]

#make column with CellID information for VR
CellID <- ID[brain_index]

#make matrix that contains the final brain points and the CellID for VR
brain <- matrix(data = c(CellID,first_brain[,"x"][brain_index],first_brain[,"y"][brain_index],first_brain[,"z"][brain_index])
                , nrow = 28461, ncol = 4, byrow = FALSE)
colnames(brain) <- c("CellID", "x", "y", "z")

#make txt file containing the reshaped matrix
write.table(brain, "braingrid", sep = "\t")

#add genes
#import selection of genes
geneselection <- read.delim("selection1.txt", header = TRUE, sep = "\t")

#make index for gene1
ind_gene <- geneselection$X733[which(geneselection$X3 == 3)]

xgene <- brain[,"x"][ind_gene]
ygene <- brain[,"y"][ind_gene]
zgene <- brain[,"z"][ind_gene]

CellIDgene <- as.character(ID[ind_gene])

#make matrix for gene
gene <- matrix(data = c(CellIDgene,xgene,ygene,zgene), nrow = length(ind_gene), ncol = 4, byrow = FALSE)
colnames(gene) <- c("CellID", "x", "y", "z")

#make one count matrix
allunifieddata_brain <- matrix(0L, nrow = nrow(brain), ncol = 1)
rownames(allunifieddata_brain) <- brain[,"CellID"]
colnames(allunifieddata_brain) <- c("gene")

matrixind <- which(is.na(match(rownames(allunifieddata_brain), gene[,"CellID"])) == FALSE)
allunifieddata_brain[,"gene"][matrixind] = 1
matrixind_na <- which(is.na(match(rownames(allunifieddata_brain), gene[,"CellID"])) == TRUE)

#add background noise
noise_generegions <- as.integer(runif(length(allunifieddata_brain[,"gene"][matrixind]), min = 1, max = 15))
allunifieddata_brain[,"gene"][matrixind] <- (allunifieddata_brain[,"gene"][matrixind] + noise_generegions)

Probability <- rbinom(length(allunifieddata_brain[,"gene"][matrixind_na]), 1, 0.01) 
noise_completebrain <- as.integer(runif(length(allunifieddata_brain[,"gene"][matrixind_na]), min = 0, max = 3)) * Probability
allunifieddata_brain[,"gene"][matrixind_na] <- (allunifieddata_brain[,"gene"][matrixind_na]+ noise_completebrain)

#plot gene inside the brain grid
open3d()
plot3d(brain[,"x"], brain[,"y"], brain[,"z"], col = ifelse(allunifieddata_brain[,"gene"] >= 1, "blue", "grey"), size = 6)

#kernel density plot
#make x,y and z coordinates for kernel density of the gene in allunifieddata (without background noise)
gene_k <- which(allunifieddata_brain[,"gene"]>=3)
xgene_k <- brain[,"x"][gene_k]
ygene_k <- brain[,"y"][gene_k]
zgene_k <- brain[,"z"][gene_k]

#compute kernel densities
kernel_brain <- kde3d(brain[,"x"], brain[,"y"], brain[,"z"], n = 100)
brain_triangles<- contour3d(kernel_brain$d, exp(-12), kernel_brain$x, kernel_brain$y, kernel_brain$z, color = "grey", alpha = 0.3, draw = FALSE)

kernel_gene <- kde3d(xgene_k, ygene_k, zgene_k, n = 100)
gene_triangles <- contour3d(kernel_gene$d, exp(-12), kernel_gene$x, kernel_gene$y, kernel_gene$z, color = "blue", alpha = 0.3, draw = FALSE)

#plot all kernel densities together
drawScene.rgl(list(brain_triangles, gene_triangles))

#make cellexal Objects
#Brain
BrainTriangles = NULL
for (i in 1:nrow(brain_triangles$v1)){
  
  BrainTriangles <- rbind(BrainTriangles, brain_triangles$v1[i,], brain_triangles$v2[i,], brain_triangles$v3[i,])
  
}

write.table(BrainTriangles, "BrainTriangles", sep = "\t") 

#make GIF of the current plot in the rgl window
Angle1 <- 5
gif.delay <- 3

Angle <- rep(Angle1 * pi / 180, 360/Angle1)

Animation.dir <- paste(getwd(), "/animations", sep="")

for (i in seq(Angle)) {
  view3d(userMatrix = rotate3d(par3d("userMatrix"),
                               Angle[i], 0, 1, 0))
  
  rgl.snapshot(filename=paste(paste(Animation.dir, "frame-", sep=""),
                              sprintf("%03d", i), ".png", sep=""))
}


