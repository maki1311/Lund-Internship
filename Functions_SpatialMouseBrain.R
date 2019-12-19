library("rgl")
library("misc3d")

#put everything into functions
#function to load in data by datapattern and store it in one matrix

loaddata <- function(datapattern){
  
  raw_data <- list()
  list_of_files <- list.files(pattern = datapattern)
  raw_data <- lapply(list_of_files, read.table)
  assign("raw_data", raw_data, .GlobalEnv)

#make one unique vector containing all the gene names 
  
  singlegenenames = NULL
  
  for (i in 1:length(raw_data)) {

    singlegenenames <- c(singlegenenames,(unique(c(names(raw_data[[i]]))))) 

  }
  
  genenames <- unique(singlegenenames)
  assign("genenames", genenames, .GlobalEnv)

    unified.list = list()
    
 
  for (i in 1:length(raw_data)){
    
    unifieddata <- matrix(0L, nrow = nrow(raw_data[[i]]), ncol = length(genenames))
    colnames(unifieddata) <- genenames
    rownames(unifieddata) <- rownames(raw_data[[i]])
    unifieddata[,colnames(raw_data[[i]])] <-  as.matrix(raw_data[[i]])
    unified.list[[i]] <-  unifieddata
        
  }
    
    assign("unified.list", unified.list, .GlobalEnv)
    
    allunifieddata <- NULL
    
    for ( i in 1:length(unified.list)) {
      
     allunifieddata <- rbind(allunifieddata, unified.list[[i]])
      
    }
   
     colnames(allunifieddata) <- genenames
    
     assign("allunifieddata", allunifieddata, .GlobalEnv)
     }

#function to split coordinates in big matrix and save in new 
splitdata <- function(data_to_split){
  
  coords.list = NULL      
  
  for (i in 1:length(data_to_split)) {
  
  x <- c("numeric", length = nrow(data_to_split[i]))
  y <- c("numeric", length = nrow(data_to_split[i]))
  
  
  for(j in 1:nrow(data_to_split[[i]])){
    
    x[j] <- as.numeric(strsplit(rownames(data_to_split[[i]]), "x")[[j]]) [1]
    y[j] <- as.numeric(strsplit(rownames(data_to_split[[i]]), "x")[[j]]) [2]
    
    coords.list[[i]] <- cbind(x, y, i)
  
    }
 
    assign("coords.list", coords.list, .GlobalEnv)
   
    }
  
  coordsbound <- NULL
  
  for (i in 1:length(coords.list)){
    
    coordsbound <- rbind(coordsbound, coords.list[[i]])
    
  }
    
    assign("coordsbound", coordsbound, .GlobalEnv)
    
  }

#function to plot the data without no expression values

plotdata <- function(name_of_gene){
  
  no_zero <- which(allunifieddata[,name_of_gene] != 0)
  
  open3d()
  Grad <- colorRampPalette(c("blue", "red"))
 
    x <- c(coordsbound[,1][no_zero]) 
    y <- c(coordsbound[,2][no_zero]) 
    z <- c(coordsbound[,3][no_zero])
    
  
plot3d(x, y, z, col = Grad(10), allunifieddata[,name_of_gene][no_zero], size = 6)  

  }

#function that combines loaddata function and splitdata function

processdata <- function(datapattern,data_to_split, name_of_gene){
  
  loaddata(datapattern)
  splitdata(data_to_split)
  plotdata(name_of_gene)
  
}

