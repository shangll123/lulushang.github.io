---
layout: default
title: Simulation
nav_order: 6
has_children: false
parent: SpatialPCA
permalink: /docs/Projects/SpatialPCA/Simulation
---

#### Prepare simulation datasets

The locations are generated from pixels in a image (LIBDsimu.jpg), which can be downloaded from [here](https://drive.google.com/drive/folders/18rwQjB3-g86A-M9xYPPJlHz60bfMABdE?usp=sharing).

```R
library("grid")
library("gridExtra")

library(jpeg)
mandrill = readJPEG("LIBDsimu.jpg", native = FALSE)
# > dim(mandrill)
# [1] 1100  984    3
# copy the image three times
mandrill.R = mandrill
mandrill.G = mandrill
mandrill.B = mandrill
# zero out the non-contributing channels for each image copy
mandrill.R[,,2:3] = 0
mandrill.G[,,1]=0
mandrill.G[,,3]=0
mandrill.B[,,1:2]=0

df = data.frame(
  red = matrix(mandrill[,,1], ncol=1),
  green = matrix(mandrill[,,2], ncol=1),
  blue = matrix(mandrill[,,3], ncol=1)
)

### compute the k-means clustering, to obtain regions from image by color
K = kmeans(df,5)

df$label = K$cluster

table(df$label)

### Replace the color of each pixel in the image with the mean 
### R,G, and B values of the cluster in which the pixel resides:
# get the coloring
colors = data.frame(
  label = 1:nrow(K$centers), 
  R = K$centers[,"red"],
  G = K$centers[,"green"],
  B = K$centers[,"blue"]
)
# merge color codes on to df
# IMPORTANT: we must maintain the original order of the df after the merge!
df$order = 1:nrow(df)
df = merge(df, colors)
df = df[order(df$order),]
df$order = NULL

# get mean color channel values for each row of the df.
R = matrix(df$R, nrow=dim(mandrill)[1])
G = matrix(df$G, nrow=dim(mandrill)[1])
B = matrix(df$B, nrow=dim(mandrill)[1])
  
# reconstitute the segmented image in the same shape as the input image
mandrill.segmented = array(dim=dim(mandrill))
mandrill.segmented[,,1] = R
mandrill.segmented[,,2] = G
mandrill.segmented[,,3] = B

RGB_label = R*B*G

# View the result
pdf("segmented.pdf")
grid.raster(mandrill.segmented)
dev.off()

#save(df, mandrill.segmented, R,G,B, file = "From_image_region.RData")

#> df[1:4,]
#       label red green blue         R         G         B
#200066     2   1     1    1 0.9997333 0.9994392 0.9992982
#200067     2   1     1    1 0.9997333 0.9994392 0.9992982
#200068     2   1     1    1 0.9997333 0.9994392 0.9992982
#200069     2   1     1    1 0.9997333 0.9994392 0.9992982

pixel_ind = c(1:dim(df)[1])
x_coor = rep(1:dim(mandrill.segmented)[2], each=dim(mandrill.segmented)[1])
y_coor = -rep(1:dim(mandrill.segmented)[1], times=dim(mandrill.segmented)[2])

data_groundtruth = data.frame(pixel_ind, x_coor,y_coor )
data_groundtruth$label = as.integer(as.factor(c(RGB_label)))

# remove background 
data_groundtruth_use = data_groundtruth[-which(data_groundtruth$label==5),]

set.seed(1234)
subsample = data_groundtruth_use[sample(1:dim(data_groundtruth_use)[1],10000,replace=F),]
# > min(subsample$y_coor)
# [1] -1591
subsample$y_coor = subsample$y_coor+1591 # make coordinates >=0
save(data_groundtruth, data_groundtruth_use,subsample, file = "simulate_spatial_cell_label.RData")

save(subsample, file = "LIBDsubsample.RData") 

```


#### Generate single cell data

Single cell data are downloaded from [GSE104276](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104276).

```R
library(Matrix)
library(readr)
library(openxlsx)

PFC_ref = read.table("GSE104276_all_pfc_2394_UMI_count_NOERCC.txt")
PFC_ref = as.matrix(PFC_ref)
celltype_ref = read.xlsx("GSE104276_readme_sample_barcode.xlsx",5)
PFC_ref_use = PFC_ref[,match(as.character(celltype_ref[,1]), colnames(PFC_ref))]
raw.data = PFC_ref_use
cell_types <- celltype_ref$cell_types
names(cell_types) <- as.character(celltype_ref[,1]) # create cell_types named list
cell_types <- as.factor(cell_types) # convert to factor data type
nUMI <- colSums(PFC_ref_use); names(nUMI) <- as.character(celltype_ref[,1]) # create nUMI named list
library(RCTD)
reference <- Reference(raw.data, cell_types, nUMI)
print(dim(reference@counts))
table(reference@cell_types)
#> table(reference@cell_types)
#       Astrocytes GABAergic neurons         Microglia           Neurons 
#               76               701                68              1057 
#              OPC        Stem cells 
#              117               290 
# > head(subsample)
#        pixel_ind x_coor y_coor label
# 858114    858114    781   1477     3
# 723450    723450    658    841     1
# 380568    380568    346    523     2
# 371366    371366    338    925     2
# 169923    169923    155   1068     1
# 466557    466557    425   1434     3

library(splatter)
cnts=as.matrix(reference@counts) 
# [1] 24153  2309
celltype=reference@cell_types
init_params <- splatEstimate(cnts[,which(reference@cell_types=="Neurons")])
save(init_params, file = "init_params_LIBD.RData")

```

#### Main functions in simulation

```R
simu = function(
  location,
  label,
  init_params,
  scenario,
  J,
  batch_facLoc,
  de_prop,
  de_facLoc,
  de_facScale,
  sim_seed,
  debug = FALSE
  )
{ 

  N <- nrow(location)
  set.seed(sim_seed)
  # 1.simulate count data
  noBatch <- ifelse(batch_facLoc == 0, TRUE, FALSE)
  params <- setParams(
    init_params,
    batchCells = rep(3 * N, 1),
    batch.rmEffect = noBatch,
    batch.facLoc = batch_facLoc,
    nGenes = J,
    group.prob = c(0.25, 0.25, 0.25,0.25),
    out.prob = 0,
    de.prob = de_prop,
    de.facLoc = de_facLoc,
    de.facScale = de_facScale,
    seed = sim_seed)

  sim_groups <- splatSimulate(
    params = params,
    method = "groups",
    verbose = FALSE)

  # remove cells having no expressed genes
    idx_zerosize <- apply(counts(sim_groups), MARGIN = 2, sum) == 0
    sim_groups <- sim_groups[, !idx_zerosize]

    if(scenario == 1){
      prop <- c(0.85, 0.05, 0.05,0.05)
    }else if(scenario == 2){
      prop <- c(0.45, 0.45, 0.05,0.05)
    }else if(scenario == 3){
      prop <- c(0.60, 0.30, 0.05,0.05)
    }else if(scenario == 4){
      prop <- c(0.35, 0.30, 0.30,0.05)
    }

  # 3.generate cell types
  print("generate cell types")

    ztrue <- label
    c_simu <- rep(NA, length(ztrue))
    if(scenario == 1){ # scenario 1: 4 cell types
      print("scenario == 1")
      c_simu <- ztrue # cell type assignment same as region assignment
    }else if(scenario != 1){
      print("scenario != 1")
      for(z in unique(label)){
        zi_idx <- ztrue == z # zi_idx is index for region z
        c_simu[zi_idx] <- sample(map_z2c(z), sum(zi_idx), prob = prop, replace = T)
      }
    }

    # 4.assign count data
    groups <- as.data.frame(colData(sim_groups)) %>% filter(Batch == "Batch" %&% 1)
    sim_cnt <- array(NA, c(J, N))
    for(c in 1:4){
      c_size <- sum(c_simu == c)  # true number of cells in cell type c
      c_cells <- groups$Cell[grepl(c, groups$Group)] # generated cells in group c
      cells_select <- sample(as.character(c_cells), c_size, replace = F)
      # sample same number of group c cells in real data from generated cells
      sim_cnt[, c_simu == c] <- as.matrix(counts(sim_groups)[, cells_select])
      # for positions of original cell type c cells, assign generated counts of group c
    }
    colnames(sim_cnt) <- "Cell" %&% 1:N
    rownames(sim_cnt) <- "Gene" %&% 1:J

  return(list(sim_cnt, c_simu, sim_seed))
}

map_z2c = function(z)
{
  case_when(
    z == 1 ~ c(1, 2, 3, 4),
    z == 2 ~ c(2, 3, 4, 1),
    z == 3 ~ c(3, 4, 1, 2),
    z == 4 ~ c(4, 1, 2, 3)
    )
}

```

#### Simulate data


##### 
```R
args <- as.numeric(commandArgs(TRUE))
i = args[1] # scenario
j = args[2] # repeat
print(i)

  library(tidyverse)
  library(Giotto)
  library(scater)
  library(Seurat)
  library(mclust)
  library(SC3)
  library(BayesSpace)
  library(gtools)
  library(splatter)
  library(reticulate)
  library(mclust )
  library(igraph)
  library(assertthat)
  library(SpatialPCA)
  library(ggplot2)
  
  
load("LIBDsubsample.RData")
load("init_params_LIBD.RData")
# These two R object can be downloaded from https://drive.google.com/drive/folders/18rwQjB3-g86A-M9xYPPJlHz60bfMABdE?usp=sharing.


res = simu(location=subsample[,2:3],label = subsample$label,init_params,
    scenario=i,J=5000, batch_facLoc=0, de_prop=0.5, de_facLoc=0.5, de_facScale=0.5,sim_seed=j, debug = FALSE)


# subspot level
count_mat = res[[1]]
truth = subsample$label
location=as.matrix(subsample[,2:3])
# truth = LIBDsimu_pseudo$info$z
location=as.matrix(location)
celltypes = res[[2]]
# first generate subspot level data
grid_subspot = make_grid(square_size = 4,location)

# square_size=4, 5077 spots generated
# square_size=4, 5077 spots generated
# square_size=4, 5077 spots generated


count_location_subspot = make_spot(grid_subspot,count_mat,celltypes,subsample$label)
count_subspot = count_location_subspot$count_spot
location_subspot = count_location_subspot$pseudo_location_spot/3
truth_subspot = count_location_subspot$subspottruth
spot_level_data = make_spot_from_subspot(count_location_subspot)
truth_subspot[spot_level_data$used_subspots]
# spot level
truth=spot_level_data$truth_spot
truth_empty = which(truth=="empty")
count_spot = spot_level_data$count_spot[,-truth_empty]
location_spot = spot_level_data$location_spot[-truth_empty,]
truth=truth[-truth_empty]
rownames(location_spot) = colnames(count_spot) = paste0("spot",1:dim(count_spot)[2])
rownames(count_spot) = paste0("gene",1:dim(count_spot)[1])

###############
# Run SpatialPCA
###############

LIBDsimu = CreateSpatialPCAObject(counts=count_spot, location=location_spot, project = "SpatialPCA",gene.type="spatial",sparkversion="spark", gene.number=3000,customGenelist=NULL,min.loctions = 20, min.features=20)
LIBDsimu = SpatialPCA_buildKernel(LIBDsimu, kerneltype="gaussian", bandwidthtype="Silverman")
LIBDsimu = SpatialPCA_EstimateLoading(LIBDsimu,fast=TRUE,SpatialPCnum=20)
LIBDsimu = SpatialPCA_SpatialPCs(LIBDsimu, fast=TRUE)

# Collect results
SpatialPCA_result = list()
SpatialPCA_result$SpatialPCs  = LIBDsimu@SpatialPCs
SpatialPCA_result$normalized_expr  = LIBDsimu@normalized_expr
SpatialPCA_result$location = LIBDsimu@location
pred_cluster= louvain_clustering(4,SpatialPCA_result$SpatialPCs,500 )
spotlist=rownames(SpatialPCA_result$location)
dist = as.matrix(dist(SpatialPCA_result$location))
pred_refine = refine_cluster_10x(pred_cluster, SpatialPCA_result$location, dist,shape="square") 
SpatialPCA_result$pred_cluster = pred_cluster
SpatialPCA_result$clusterlabel = pred_refine
SpatialPCA_result$truth = truth[match(rownames(LIBDsimu@location),rownames(location_spot))]
SpatialPCA_result$ARI = adjustedRandIndex(SpatialPCA_result$clusterlabel,SpatialPCA_result$truth)
SpatialPCA_result$NMI = compare(as.factor(SpatialPCA_result$clusterlabel),as.factor(SpatialPCA_result$truth), method = "nmi")
SpatialPCA_result$CHAOS = fx_CHAOS(SpatialPCA_result$clusterlabel, SpatialPCA_result$location)
SpatialPCA_result$PAS = fx_PAS(SpatialPCA_result$clusterlabel, SpatialPCA_result$location)

save(SpatialPCA_result, file = paste0("spotlevel_SpatialPCA_spatialgene_result_scenario_",i,"_rep_",j,".RData"))

###############
# BayesSpace
###############

# filter out spots with 0 counts
ind_keep=which(colSums(count_spot) > 0)
location_spot_bayesSpace=location_spot[ind_keep,]
count_spot_BayesSpace=count_spot[,ind_keep]
colnames(location_spot_bayesSpace) <- c("row", "col")

sce_LIBDsimu_ST <- SingleCellExperiment(assays = list(counts = count_spot_BayesSpace), colData = location_spot_bayesSpace)
sce_LIBDsimu_ST <- spatialPreprocess(sce_LIBDsimu_ST, platform="ST",n.PCs = 15, n.HVGs = 2000, log.normalize = T)
sce_LIBDsimu_ST <- spatialCluster(sce_LIBDsimu_ST, q=4, d=15, platform='ST',nrep=10000, gamma=3, save.chain=FALSE) 
sce_labels=sce_LIBDsimu_ST$spatial.cluster

BayesSpace_ST_result = list()
BayesSpace_ST_result$sce_LIBDsimu_ST = sce_LIBDsimu_ST
BayesSpace_ST_result$clusterlabel = sce_labels
BayesSpace_ST_result$location = location_spot_bayesSpace
BayesSpace_ST_result$truth = truth[ind_keep]
BayesSpace_ST_result$ARI = adjustedRandIndex(BayesSpace_ST_result$clusterlabel,BayesSpace_ST_result$truth)
BayesSpace_ST_result$NMI = compare(BayesSpace_ST_result$clusterlabel,BayesSpace_ST_result$truth, method = "nmi")
BayesSpace_ST_result$CHAOS = fx_CHAOS(BayesSpace_ST_result$clusterlabel, BayesSpace_ST_result$location)
BayesSpace_ST_result$PAS = fx_PAS(BayesSpace_ST_result$clusterlabel, BayesSpace_ST_result$location)

print(BayesSpace_ST_result$ARI)

save(BayesSpace_ST_result,  file = paste0("spotlevel_BayesSpace_hvggene_result_scenario_",i,"_rep_",j,".RData"))


###############
# prepare SpaGCN data
# and then run in python
###############

print("SpaGCN")
# SpaGCN uses all genes

library(DropletUtils)
write10xCounts(
  paste0("Spotlevel_SpaGCN_scenairo_",i,"_rep_",j,"_count_allgene.h5"),
  count_spot,
  barcodes = colnames(count_spot),
  gene.id = rownames(count_spot),
  gene.symbol = rownames(count_spot),
  gene.type = "Gene Expression",
  overwrite = FALSE,
  type = c( "HDF5"),
  genome = "unknown",
  #version = c("2", "3"),
  #chemistry = "Single Cell 3' v3",
  original.gem.groups = 1L,
  library.ids = "custom"
)


```










