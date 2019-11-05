# Script to identifying edge clusters in APT datasets
# Required input of POS File, Cluster Range F
# Find clusters on the edge of an alpha hull
#   Input: pos file and cluster stats file path
#   Cluster stats may be exported from IVAS (csv)
#   AtomicDensity = Atomic density of material (atoms/nm3)
#   DetectionEfficiency = Detection efficiency of leap used
#   SamplingFraction = How much one would like to sample the POS file by
#   Returns a list of which cluster IDs are deemed to be on the edge of the dataile, Cluster Analysis CSV
# Written by Ben Jenkins and Andy London

setwd(dirname(parent.frame(2)$ofile))
require("tidyverse")
require("geometry")
require("alphashape3d")
require("spatstat")

#These functions allow one to import a pos file (sampled or not) into R with x, y, z, and mass-to-charge-ratio information
source("Test Files\\read.pos.sampled.R")
source("Test Files\\read.pos.R")
source("Test Files/writeposR.R")
#### Function for identifying edge clusters ####
findEdgeClustersConvex <- function (posFileName, 
                                    clusterStatsFile,
                                    clusterIndexPosFile, 
                                    clusterPosFile,
                                    SamplingFraction,
                                    NNDMultiplier) {
  
  
  
  # get time for timing execution
  StartTime <- Sys.time()
  
  #### Read the cluster statistics file ####
  #cluster file path
  ClusterImport <- read.csv(clusterStatsFile, skip = 10)
  
  #### Cluster Import ####
  ClusterImport <-
    ClusterImport %>% select(
      id=X,
      Solute.Ions,
      Ranged.Ions,
      Total.Ions,
      Center_x..nm..Ranged,
      Center_y..nm..Ranged,
      Center_z..nm..Ranged
    ) %>%
    mutate(
      xpos = Center_x..nm..Ranged,
      ypos = Center_y..nm..Ranged,
      zpos = Center_z..nm..Ranged
    ) %>%
    filter(grepl("Cluster", id))
  
  ClusterImportedData <- ClusterImport
  
  # load index cluster file to define cluster centres and boundary by convex hull
  clrIndxData <- read.pos(clusterIndexPosFile)
  ClusterPosFile <- read.pos(clusterPosFile)
  nclustersDet <- max(clrIndxData[,4])
  clrHull_list <- vector("list", nclustersDet)

  for (i in 1:nclustersDet) {
    thisClr <- clrIndxData[clrIndxData[,4]==i,1:3]
    COM <- summarise_all(thisClr,"mean")
    tempList<-convhulln(thisClr)
    clrHull <- thisClr[unique(array(tempList)),1:3]
    clrHull$id <- i
    clrHull_list[[i]] <- clrHull
  }
  clrHull_all <- bind_rows(clrHull_list)
  
  
  #### Load filtered pos file to improve speed ####
  FilterPosFile <- read.pos.sampled(posFileName, SamplingFraction)
  
  #### Parameters For Calculating Alpha Value####
  AlphaValue <<- NNDMultiplier*round(ceiling(100*(max(nndist(FilterPosFile %>% select(x,y,z), k=1)))),2)/100
  
  print(paste0(
    "Alpha Value: ",
    AlphaValue,
    " Sampling fraction: ",
    SamplingFraction
  ))
  
  #### Using alpha-shape 3d ####
  M <- as.matrix(FilterPosFile[, 1:3])
  
  # free some memory
  rm(FilterPosFile)
  gc()
  
  AlphaHullShape <<- ashape3d(M, alpha =  c(AlphaValue, 10))

  #### Ensuring there is only one connected volume ####
  comp <- components_ashape3d(AlphaHullShape, indexAlpha = "all")
  NumberVolumes <- as.data.frame(table(comp[[1]]))
  
  
  #### Determing which clusters are edge clusters and returning values ####
  if (nrow(NumberVolumes) != 1) {
    stop(paste0("Error. Multiple volumes created"))
  } else{
    print(paste0("Alpha value has created one volume"))
    
    # determining which clusters are edge clusters
    ClustersInHull <-
      inashape3d(AlphaHullShape, indexAlpha = 1, as.matrix(clrHull_all[,1:3]))
    
    # getting names of clusters that are edge
    ClusterLocation2 <- cbind(clrHull_all, ClustersInHull) %>% 
      dplyr::filter(ClustersInHull==FALSE) %>% 
      group_by(id) %>% 
      summarise(ID=first(id))
    
    TotalEdgeClusters <- sort(as.numeric(as.matrix(ClusterLocation2[,1])))
  }
  
  
  print(paste0(
    "Toal number of clusters detected: ",
    length(TotalEdgeClusters)
  ))
  NumberOfAssignedEdgeClusters <<- length(TotalEdgeClusters)
  # Pritn total exe time
  EndTime <- Sys.time()
  TotalTime <<- EndTime - StartTime
  print(TotalTime)
  
  EdgeClusters <<- TotalEdgeClusters
  NumberClusters <<- max(as.numeric(gsub("Cluster ","",ClusterImportedData$id)))
  #return(list(TotalEdgeClusters=TotalEdgeClusters, ClusterImportedData=ClusterImportedData))
  
  # Return pos file containing only ions in edge clusters and only ions not in edge clusters
  IonsInEdgeClusters <- vector()
  for (i in EdgeClusters) {
    thisClr <- clrHull_all[clrHull_all[,4]==i,1:3]
    COM <- summarise_all(thisClr,"mean")
    tempList<-convhulln(thisClr)
    IonsInEdgeClustersi <- which(inhulln(tempList, as.matrix(ClusterPosFile %>% select(-m))))
    IonsInEdgeClusters <- append(IonsInEdgeClustersi,IonsInEdgeClusters)
  }
  
  EdgeClusterPos <- ClusterPosFile %>% slice(IonsInEdgeClusters)
  NonEdgeClusterPos <- setdiff(ClusterPosFile,EdgeClusterPos)
  writeposR(EdgeClusterPos, "Edge Clusters Ranged.pos")
  writeposR(setdiff(ClusterPosFile,EdgeClusterPos), "Non Edge Clusters Ranged.pos")
}

#### Parameters For Calculating Alpha Value####
# How much one would like to sample the POS file by
SamplingFraction <- 0.005 
# What NND multiplier should be used to calculate alpha 
AlphaNND <- 4

#### File paths for fundtion ####
findEdgeClustersConvex("Test Files/006_full.pos",
                       "Test Files/006_full - Top-Level ROI - Cluster Analysis (Ni, Cu).csv",
                       "Test Files/006_full.cluster.indexed.pos",
                       "Test Files/006_full.cluster.pos",
                       SamplingFraction,
                       AlphaNND)

# View Alpha shape that is generated
plot(AlphaHullShape, indexAlpha = 1)

# Calculate number density (assuming correction factor of 0.5 and that alpha volume generated around point cloud is correct)
NumberDensity = (NumberClusters - (0.5 * NumberOfAssignedEdgeClusters))/
  (volume_ashape3d(AlphaHullShape, byComponents = FALSE, indexAlpha = 1)*((10^-9)^3))
print(paste0("The number density of clusters in this dataset is: ", NumberDensity,
             ". A correction factor of 0.5 has been used and the following assumption has been made: the volume generated by alpha shape is the correct volume of your dataset."))
