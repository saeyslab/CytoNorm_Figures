## module load gridengine
## qrsh -now no -l h_vmem=10G,h_rt=1:0:0
## module unload gcc
## export LD_LIBRARY_PATH="/software/shared/apps/x86_64/gcc/4.8.0/lib/:/software/shared/apps/x86_64/gcc/4.8.0/lib64:$LD_LIBRARY_PATH"
## module load R/x86_64/3.5.1
## R
##   

taskID <- as.numeric(commandArgs(trailingOnly=TRUE)[1])-1
param1 <- (taskID %% 5) + 1
param2 <- (floor(taskID / 5)) + 1

####################
## Load libraries ##
####################

library(flowCore)
library(FlowSOM)
library(CytoNorm)

#########################
## Set some parameters ##
#########################

dir <- "/group/irc/personal/sofievg/Stanford/OriginalData/GatesExtraControls"
results_dir <- "/group/irc/personal/sofievg/Stanford/NormalizedData/GatesExtraControls"
norm_dir <- file.path(results_dir, taskID)
dir.create(norm_dir)

setwd(norm_dir)

files <- list.files(dir, pattern="Unstim_Control.*.fcs$")

train_files <- grep("Unstim_Control_1", files, value = TRUE)
train_labels <- gsub("_[UI].*", "", gsub("Gates_", "", train_files))

test_files <- grep("Unstim_Control_2", files, value = TRUE)
test_labels <- gsub("_[UI].*", "", gsub("Gates_", "", test_files))

norm_markers <- c(48, 46, 43, 45, 20, 16, 21, 19, 22, 50, 47, 40, 44, 33,  # 17,
                 11, 18, 51, 14, 23, 32, 10, 49, 27, 24, 31, 42, 37, 39, 34, 
                 41, 26, 30, 28, 29, 25, 35)

# Use description for marker name when available
ff <- read.FCS(file.path(dir, train_files[1]))
marker_names <- get_markers(ff, colnames(ff))

surface_markers <- c(grep("CD", marker_names),
                     grep("CC", marker_names),
                     grep("HLADR", marker_names),
                     grep("TCRgd", marker_names))
surface_markers <- surface_markers[which(surface_markers %in% norm_markers)]

# Build transformlist for further use
transform_list <- transformList(colnames(ff)[norm_markers],
                                cytofTransform)

transform_list.reverse <- transformList(colnames(ff)[norm_markers], 
                                        cytofTransform.reverse)

############################
## Load the manual gating ##
############################

manual <- readRDS("/group/irc/personal/sofievg/Stanford/ManualGating/manual_GatesExtraControls.RDS")
cell_types <- c("All", levels(manual[[1]]$manual)[-1])

###################
## Normalization ##
###################

#for(nQ in c(51,21,6,11,101)){
#  gridSize <- c(15,10,20,5,25)[param1]
#  #for(gridSize in c(15,10,20,5,25)){
#  nClus <- c(25,15,35,5,45)[param2]
#  #for(nClus in c(25,15,35,5,45)){

for(nClus in c(55,65,75)){
  gridSize <- c(15,10,20,5,25)[param1]
  nQ <- c(51,21,6,11,101)[param2]
   if(nClus < gridSize*gridSize){
    setting <- paste0("nClus_",nClus,"_gridSize_",gridSize,"_nQ_",nQ)
    print(setting)
    
    for(i in 1:3){
      # Learn the splines
      t_train <- system.time(capture.output(
        normModel <- CytoNorm.train(files = file.path(dir, train_files),
                                     labels = train_labels,
                                     channels = names(marker_names)[norm_markers],
                                     transformList = transform_list,
                                     FlowSOM.params = list(nCells = 1000000,
                                                           nClus = nClus,
                                                           xdim = gridSize,
                                                           ydim = gridSize),
                                     normMethod.train = QuantileNorm.train,
                                     normParams = list(nQ = nQ),
                                     seed = i)
      ))
      # Normalize the data
      t_norm <- system.time(capture.output(
        CytoNorm.normalize(normModel,
                            files= file.path(dir, test_files),
                            labels = test_labels,
                            transformList = transform_list,
                            transformList.reverse = transform_list.reverse,
                            outputDir = norm_dir)
      ))
      
      
      evaluation_norm <- emdEvaluation(file.path(norm_dir, 
                                                 paste0("Norm_",test_files)),
                                       transformList = transform_list,
                                       manual = lapply(manual, function(x){x$manual}),
                                       markers = names(marker_names)[norm_markers])
      
      save(evaluation_norm, t_train, t_norm,
           file = file.path(results_dir,paste0("eval",i,"_",setting,".Rdata")))
    }
  }
  #}
  #}
}