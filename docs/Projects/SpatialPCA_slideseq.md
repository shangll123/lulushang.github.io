---
layout: default
title: Slideseq analysis
nav_order: 4
has_children: false
parent: SpatialPCA
permalink: /docs/Projects/SpatialPCA/Slideseq
---


##### We illustrate the benefits of SpatialPCA through four different downstream analyses: spatial transcriptomics visualization, spatial domain detection, spatial trajectory inference on the tissue, and high-resolution spatial map reconstruction. 

### Section 1: Data processing.
Load packages.
```R
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(plotly)
  library(STutility)
  library(zeallot)
  library(openxlsx)
```

We used SPARK-X to select spatial genes as input for SpatialPCA. The "Puck_180430_6_zero_removed_nomt.rds" data could be downloaded from [here](https://drive.google.com/file/d/1t1ITBO6RvlV_GK5ze6RYC59lzRIf0n3L/view?usp=sharing). 
```R
# Calling SPARK-X functions. 
source("~/kdc_batch_eigen_correct.R")
sourceCpp("~/kdc.cpp")
# For the details of the SPARK-X method please refer to https://xzhoulab.github.io/SPARK.

load(paste0("Puck_180430_6_zero_removed_nomt.rds"))
slideseqSeu1 <- CreateSeuratObject(counts = sp_count, project = "slideseq", min.cells = 20, min.features = 20)
slideseqSeu1 = SCTransform(slideseqSeu1, 
              return.only.var.genes = FALSE, 
              variable.features.n = NULL, 
              variable.features.rv.th = 1.3)
sp_count1 = sp_count[match(rownames(slideseqSeu1@assays$SCT@scale.data),rownames(sp_count)),match(colnames(slideseqSeu1@assays$SCT@scale.data),colnames(sp_count))]
raw_loc1 = location[match(colnames(sp_count1), rownames(location)),]
t1 <- proc.time()
KDC1 <- kdc_mk_spark(sp_count1,raw_loc1,numCores=10,filter=F,option="mixture",verbose=FALSE)
t2 <- proc.time() - t1
runTime_kdc = t2

# sum(KDC1$res_mtest$adjustedPval<0.1)
# sum(KDC1$res_mtest$adjustedPval<0.05)
# sum(KDC1$res_mtest$adjustedPval<0.01)

# > sum(KDC1$res_mtest$adjustedPval<0.1)
# [1] 914
#> sum(KDC1$res_mtest$adjustedPval<0.05)
#[1] 787
# > sum(KDC1$res_mtest$adjustedPval<0.01)
# [1] 522

expr = slideseqSeu1@assays$SCT@scale.data[which(KDC1$res_mtest$adjustedPval<0.05),]
info = scale(raw_loc1)

```



### Section 2: Run SpatialPCA.

SpatialPCA:
```R
expr=scale_expr(expr)
dat = data_prepare_func(expr, info)
bandwidth = bandwidth_select(expr, info,method="Silverman")
K=kernel_build(kernelpara="gaussian", dat$ED2,bandwidth) 

# Set the number of PCs.
require(RSpectra)
PC_num = 20
dat$YMt = t(dat$YM)
dat$KYM = K%*%dat$YMt
dat$K = K
eigen_res = eigs_sym(K, k=PC_num, which = "LM")
dat$delta = eigen_res$values
dat$U = eigen_res$vectors

# In this slideseq data we have ~20k locations. 
# Below first step takes about 10min, 10Gb
Est_para = SpatialPCA_estimate_parameter_largedata(dat_input=dat,PCnum=PC_num) 
# Below second step takes about 25s, 3Gb
Est_W = SpatialPCA_estimate_W_largedata(Est_para$par, dat,PCnum=PC_num)
# Below third step takes about 3h, 17Gb, mainly spent on taking an inverse of a 20k by 20k matrix.
Est_Z = SpatialPCA_estimate_Z_largedata(Est_para$par,dat,Est_W,PCnum=PC_num)
```

PCA and NMF:
```R
dat$Z_pca=get_PCA(expr, PCnum=20)
dat$Z_NMF = get_NMF(expr, PCnum=20)
```

HMRF:
```python
# prepare files in r
setwd("~/slideseq/HMRF")
ids = c(1:dim(expr)[2])
exprdata = as.data.frame(t(expr))
writeexpr = cbind(ids,exprdata)
write.table(writeexpr, file = "~/slideseq/HMRF/expr_marker.txt", row.names=F,col.names=F, quote=F )
writelocation = cbind(ids,0,info )
write.table(writelocation, file = "~/slideseq/HMRF/location.txt", row.names=F,col.names=F, quote=F )
gene = rownames(expr)
write.table(gene, file = "~/slideseq/HMRF/genes", row.names=F,col.names=F, quote=F )

# run HMRF in python
import timeit
import sys
import math
import os
import numpy as np
import pandas as pd
import smfishHmrf.reader as reader
from smfishHmrf.HMRFInstance import HMRFInstance
from smfishHmrf.DatasetMatrix import DatasetMatrix, DatasetMatrixSingleField
import smfishHmrf.visualize as visualize
import smfishHmrf.spatial as spatial

directory = "~/slideseq/HMRF"
# read input files
genes=reader.read_genes("%s/genes" % directory)
Xcen,field=reader.read_coordinates("%s/location.txt" % directory)
expr=reader.read_expression_matrix("%s/expr_marker.txt" % directory)
# create a DatasetMatrixSingleField instance with the input files
this_dset = DatasetMatrixSingleField(expr, genes, None, Xcen)
# compute neighbor graph (first test various distance cutoff: 0.3%, 0.5%, 1% of lowest pairwise distances)
this_dset.test_adjacency_list([0.3, 0.5, 1], metric="euclidean")
this_dset.calc_neighbor_graph(0.3, metric="euclidean")
this_dset.calc_independent_region()
new_genes = reader.read_genes("%s/genes" % directory)
new_dset = this_dset.subset_genes(new_genes)
os.getcwd() 
outdir = "spatial.HMRF.cluster"
os.mkdir(outdir)
start = timeit.default_timer()
os.chdir('~/slideseq/HMRF') 
this_hmrf = HMRFInstance("slideseq", outdir, new_dset, 8,  (0, 0.5, 30), tolerance=1e-10)
this_hmrf.init(nstart=1000, seed=-1)
this_hmrf.run()
stop = timeit.default_timer()
print('Time: ', stop - start)  
visualize.domain(this_hmrf, 8, 10, dot_size=8, size_factor=10, outfile="cluster8.visualize.beta.%.1f.png" % 5.0)


# collect results
k = 8
name = "slideseq"
betas = seq(from=0,to=14.5,by=0.5)
python_path = "/net/mulan/home/shanglu/anaconda3/bin/python3.7"
get_result_path = system.file("python", "get_result2.py", package = 'Giotto')
output_data = "~/her2st/HMRF/spatial.HMRF.cluster"
cluster_HMRF = list()
for(num in 1:length(betas)){
b = betas[num]
print(num)
    result_command = paste0(python_path, ' ', get_result_path,
                            ' -r ', output_data,
                            ' -a ', name,
                            ' -k ', k,
                            ' -b ', b)
    print(result_command)
    output = system(command = result_command, intern = T)
    cluster_HMRF[[num]] = data.table::data.table(temp_name = output)$temp_name

}

```


### Section 3: Spatial transcriptomics visualization.

Visualization by PC values.
```R
p_PCs = plot_factor_value(info,dat$Z_spatial,textmethod="SpatialPCA",pointsize=0.5,textsize=15)
ggarrange(p_PCs[[1]], p_PCs[[2]], p_PCs[[3]],
          ncol = 3, nrow = 1)

p_PCs = plot_factor_value(info,dat$Z_pca,textmethod="PCA",pointsize=0.5,textsize=15)
ggarrange(p_PCs[[1]], p_PCs[[2]], p_PCs[[3]],
          ncol = 3, nrow = 1)

p_PCs = plot_factor_value(info,dat$Z_NMF,textmethod="NMF",pointsize=0.5,textsize=15)
ggarrange(p_PCs[[1]], p_PCs[[2]], p_PCs[[3]],
          ncol = 3, nrow = 1)

```
Visualization by RGB plots.
```R
library(ggpubr)

p1 = plot_RGB_tSNE(info,dat$Z_spatial,pointsize=0.5)
p2 = plot_RGB_tSNE(info,dat$Z_pca,pointsize=0.5)
p3 = plot_RGB_tSNE(info,dat$Z_NMF,pointsize=0.5)
pdf("RGB_tSNE.pdf",width=15, height=5)
ggarrange(p1, p2, p3, 
           labels = c("SpatialPCA", "PCA", "NMF"),
          ncol = 3, nrow = 1)
dev.off()

p1 = plot_RGB_UMAP(info,dat$Z_spatial,pointsize=0.5)
p2 = plot_RGB_UMAP(info,dat$Z_pca,pointsize=0.5)
p3 = plot_RGB_UMAP(info,dat$Z_NMF,pointsize=0.5)
pdf("RGB_UMAP.pdf",width=15, height=5)
ggarrange(p1, p2, p3, 
          labels = c("SpatialPCA", "PCA", "NMF"),
          ncol = 3, nrow = 1)
dev.off()

```

### Section 4: Use clustering results from SpatialPCA to obtain spatial domains.
We used Louvain clustering results in the main analysis to save time, because Walktrap in large data is pretty slow. In slideseq data we used 8 clusters. 
```R
N_k_nearest=c(seq(from=30,to=90,by=5),100,150,200,250)
louvain_cluster_SpatialPCA = louvain_clustering(knearest = N_k_nearest, dat$Z_spatial)
N_k_nearest=c(seq(from=100,to=1000,by=50))
louvain_cluster_PCA = louvain_clustering(knearest = N_k_nearest, dat$Z_pca)
N_k_nearest=c(seq(from=5000,to=10000,by=500))
louvain_cluster_NMF = louvain_clustering(knearest = N_k_nearest, dat$Z_NMF)

colnames(info) = c("sdimx","sdimy")
meta_data = info
meta_data$SpatialPCA_Louvain = louvain_cluster_SpatialPCA$cluster_label[[28]]
meta_data$PCA_Louvain = louvain_cluster_PCA$cluster_label[[5]]
meta_data$NMF_Louvain = louvain_cluster_NMF$cluster_label[[11]]
meta_data$HMRF = cluster_HMRF[[21]]

```

### Section 5: Clustering results visualization.

```R
loc1 = info[,1]
loc2 = info[,2]
color_use_SpatialPCA = c("#66C2A5", "lightyellow2", "cornflowerblue" ,"#E78AC3", "skyblue1" ,"#FFD92F" ,"lightcyan2", "coral")
pdf("Plot_8_clusters_SpatialPCA_Louvain.pdf",width=5,height=5)
plot_cluster(loc1,loc2,meta_data$SpatialPCA_Louvain,pointsize=0.7,text_size=10 ,"",color_use_SpatialPCA)
dev.off()

color_use_PCA = c("#66C2A5",  "cornflowerblue" ,"#FFD92F","#E78AC3","lightcyan2","lightyellow2", "skyblue1" ,  "coral")
pdf("Plot_8_clusters_PCA_Louvain.pdf",width=5,height=5)
plot_cluster(loc1,loc2,meta_data$PCA_Louvain,pointsize=0.7,text_size=10 ,"",color_use_PCA)
dev.off()

color_use_NMF = c("#66C2A5", "lightyellow2", "cornflowerblue" ,"#E78AC3", "skyblue1" ,"#FFD92F" ,"lightcyan2", "coral")
pdf("Plot_8_clusters_NMF_Louvain.pdf",width=5,height=5)
plot_cluster(loc1,loc2,meta_data$NMF_Louvain,pointsize=0.7,text_size=10 ,"",color_use_NMF)
dev.off()

color_use_HMRF = c("#66C2A5", "lightyellow2", "cornflowerblue" ,"#E78AC3", "skyblue1" ,"#FFD92F" ,"lightcyan2", "coral")
pdf("Plot_8_clusters_HMRF.pdf",width=5,height=5)
plot_cluster(loc1,loc2,meta_data$HMRF,pointsize=0.7,text_size=10 ,"",color_use_HMRF)
dev.off()


```
### Section 6: Trajectory analysis.

Focus on granule cell layer:
```R
library(slingshot)
library(SingleCellExperiment)


cellid = which(meta_data$SpatialPCA_Louvain %in% c("1","6","3"))
sim_spatialPCA_granule = SingleCellExperiment(assays = list(counts = result$Z_spatial[,cellid]))
reducedDims(sim_spatialPCA_granule) <- SimpleList(DRM = t(result$Z_spatial[,cellid]))
colData(sim_spatialPCA_granule)$louvain <- factor(meta_data$SpatialPCA_Louvain[cellid])    
sim_spatialPCA_granule  <-slingshot(sim_spatialPCA_granule, clusterLabels = 'louvain', reducedDim = 'DRM',start.clus="1" )
summary(sim_spatialPCA_granule@colData@listData)
# > summary(sim_spatialPCA_granule@colData@listData)
#                   Length Class  Mode   
# louvain           8997   factor numeric
# slingPseudotime_1 8997   -none- numeric
save(sim_spatialPCA_granule, file = "slideseq_slingshot_sim_SpatialPCA_granule.RData")
spatialPCA_granule_cluster = meta_data$SpatialPCA_Louvain[cellid]


method = "SpatialPCA"
gridnum = 25
arrowsize = 0.7
arrowlength = 0.2
pointsize=2

sim=sim_spatialPCA_granule
color_in = color_use_SpatialPCA[c(1,3,6)]
color_in[4] = "black"
clusterlabels = spatialPCA_granule_cluster

for(traj in 1:(length(sim@colData@listData)-1)){
pseudotime_use = sim@colData@listData[[traj+1]]
p_list = plot_trajectory(pseudotime_use, info[cellid,],clusterlabels,gridnum,color_in,pointsize=pointsize ,arrowlength=arrowlength,arrowsize=arrowsize,textsize=30)
pdf(paste0("Granule_",method,"_slingPseudotime_",traj,"_grid",gridnum,"_arrowsize",arrowsize,".pdf"))
print(p_list[[1]])
print(p_list[[2]])
print(p_list[[3]])
print(p_list[[4]])
print(p_list[[5]])
dev.off()
}


```
Trajectory inference on whole tissue section:
```R
library(slingshot)
library(SingleCellExperiment)

# on whole tissue section
method = "SpatialPCA"
sim<- SingleCellExperiment(assays = list(counts = as.matrix(result$Z_spatial)))
reducedDims(sim) <- SimpleList(DRM = t(result$Z_spatial))
colData(sim)$Louvain <- factor(meta_data$SpatialPCA_Louvain)    
sim  <-slingshot(sim, clusterLabels = 'Louvain', reducedDim = 'DRM',start.clus="7" )
summary(sim@colData@listData)
sim_SpatialPCA=sim
save(sim_SpatialPCA, file = "slideseq_sim_SpatialPCA.RData")


#> summary(sim@colData@listData)
#                  Length Class  Mode
#Louvain           20982  factor numeric
#slingPseudotime_1 20982  -none- numeric
#slingPseudotime_2 20982  -none- numeric
#slingPseudotime_3 20982  -none- numeric

pseudotime_traj1 = sim@colData@listData$slingPseudotime_1
pseudotime_traj2 = sim@colData@listData$slingPseudotime_2
pseudotime_traj3 = sim@colData@listData$slingPseudotime_3
clusterlabels = meta_data$SpatialPCA_Walktrap
gridnum = 50
color_in = c(  "plum1", "dodgerblue","mediumaquamarine",  "palegreen4", "chocolate1","lightblue2","#F0E442",  "black","#CC79A7","mediumpurple","seagreen1")
gridnum = 50
arrowsize = 1
arrowlength = 0.2
pointsize=2

for(traj in 1:(length(sim@colData@listData)-1)){
pseudotime_use = sim@colData@listData[[traj+1]]
p_list = plot_trajectory(pseudotime_use, info,clusterlabels,gridnum,color_in,pointsize=pointsize ,arrowlength=arrowlength,arrowsize=arrowsize,textsize=22)
pdf(paste0("slideseq_",method,"_slingPseudotime_",traj,"_grid",gridnum,"_arrowsize",arrowsize,".pdf"),width=10,height=10)
print(p_list[[1]])
print(p_list[[2]])
print(p_list[[3]])
print(p_list[[4]])
print(p_list[[5]])
dev.off()
}

```


### Section 7: Cell type deconvolution.
We used RCTD to perform cell type deconvolution. The reference data is from cerebellum data in [http://dropviz.org](http://dropviz.org)
```R
library(RCTD)
library(Matrix)
# make reference data
# try DropViz
# install.packages('DropSeq.util_2.0.tar.gz', repos=NULL)
library(DropSeq.util)
 dge.path <- "~/RCTD/F_GRCm38.81.P60Cerebellum_ALT.raw.dge.txt"
 dge <- loadSparseDge(dge.path) 
cell_cluster_outcomes = readRDS("~/RCTD/F_GRCm38.81.P60Cerebellum_ALT.cell_cluster_outcomes.RDS")
subcluster = readRDS("~/RCTD/F_GRCm38.81.P60Cerebellum_ALT.subcluster.assign.RDS")
clusterassign = readRDS("~/RCTD/F_GRCm38.81.P60Cerebellum_ALT.cluster.assign.RDS")
annotationBraincellAtlas = readRDS("~/RCTD/annotation.BrainCellAtlas_Saunders_version_2018.04.01.RDS")

# Cerebellum  Neuron  GranularNeuron_Gabra6 [#1]
# Cerebellum  Endothelial Endothelial_Flt1 [#10]
# Cerebellum  Fibroblast-Like Fibroblast-Like_Dcn [#11]
# Cerebellum  Neuron  PurkinjeNeuron_Pcp2 [#2]
# Cerebellum  Neuron  Interneurons_Pvalb [#3]
# Cerebellum  Neuron  Interneurons_and_Other_Nnat [#4]
# Cerebellum  Microglia_Macrophage  Microglia_Macrophage_C1qb [#5]
# Cerebellum  Oligodendrocyte_Polydendrocyte  Oligodendrocyte_Polydendrocyte_Tfr_Tnr [#6]
# Cerebellum  Astrocyte BergmannGlia_Gpr37l1 [#7]
# Cerebellum  Astrocyte Astrocyte_Gja1 [#8]
# Cerebellum  Choroid_Plexus  Choroid_Plexus_Ttr [#9]

reference_clusters = as.character(paste0(1:11))
reference_clusters_common_names = c("GranularNeuron","PurkinjeNeuron","Interneurons","Interneurons_and_Other_Nnat","Microglia_Macrophage",
  "Oligodendrocyte_Polydendrocyte","BergmannGlia","Astrocyte","Choroid_Plexus","Endothelial","Fibroblast")
datt = data.frame(reference_clusters,reference_clusters_common_names )
cell_cluster_outcomes$reason = NULL
cell_cluster_outcomes = na.omit(cell_cluster_outcomes)
raw.data = dge[,match(rownames(cell_cluster_outcomes),colnames(dge))]
cell_cluster_outcomes$nUMI = colSums(raw.data)
cell_cluster_outcomes$liger_ident_coarse = cell_cluster_outcomes$cluster
reference = Seurat::CreateSeuratObject(raw.data, meta.data = cell_cluster_outcomes)
ref_celltype = cell_cluster_outcomes$cluster
names(ref_celltype) = rownames(cell_cluster_outcomes)
reference <- Reference(raw.data, cell_types=ref_celltype, nUMI=colSums(raw.data))
load(paste0("~/Puck_180430_6_zero_removed_nomt.rds"))
SCTcount = sp_count[,match(colnames(expr),colnames(sp_count))]
location = location[match(colnames(SCTcount),rownames(location)),]
colnames(location) = c("x","y")
barcodes = rownames(location)
coords = data.frame(barcodes, "x" = location[,1], "y" = location[,2])
rownames(coords) = NULL
coords[,2:3] = scale(coords[,2:3])
write.csv(location, "BeadLocationsForR.csv",quote=F)
colnames(coords)[2] = "x"
colnames(coords)[3] = "y"
SCTcount_dataframe = as.data.frame(SCTcount)
rownames(SCTcount_dataframe) = NULL
coords = tibble::column_to_rownames(coords, var = "barcodes")
coords$barcodes <- NULL
nUMI <- colSums(SCTcount)
puck = SpatialRNA(coords, SCTcount,nUMI)

library(doParallel)
myRCTD <- create.RCTD(puck, reference, max_cores = 5)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')

results <- myRCTD@results
dim(myRCTD@results$results_df)
#> dim(myRCTD@results$results_df)
#[1] 9158    9

metadata_RCTD = myRCTD@results$results_df
meta_data$cellid = rownames(meta_data)
metadata_RCTD$cellid = rownames(meta_data_RCTD)
metadata_RCTD = merge(meta_data,metadata_RCTD, by="cellid" )
metadata_RCTD$celltype = datt$reference_clusters_common_names[match(metadata_RCTD$first_type, datt$reference_clusters)]
metadata_RCTD$celltype = as.factor(metadata_RCTD$celltype)

save(meta_data_RCTD, file = "meta_data_RCTD.RData")
save(meta_data, file = "meta_data.RData")

# > table(metadata_RCTD$celltype)

#                      Astrocyte                   BergmannGlia 
#                            178                            704 
#                 Choroid_Plexus                    Endothelial 
#                             85                             21 
#                     Fibroblast                 GranularNeuron 
#                             35                           5129 
#                   Interneurons    Interneurons_and_Other_Nnat 
#                            183                            539 
#           Microglia_Macrophage Oligodendrocyte_Polydendrocyte 
#                             43                            920 
#                 PurkinjeNeuron 
#                           1321 
                          
clusternum = 8
celltypes = 11

method = "SpatialPCA"
cbp2=c(
"#66C2A5", "#FC8D62","#8DA0CB" ,"#E78AC3", "#A6D854", 
 "#E5C494" ,"skyblue2","#1F78B4", 
"#B2DF8A", "#33A02C" ,"#FB9A99", "#E31A1C", "#BF5B17" ,
"#7FC97F", "#BEAED4" , "#FFFF99" ,"#F0027F" )

p = plot_celltype_barplot_total100(clusternum, celltypes, metadata_RCTD,method,cbp2,textsize=30)
pdf(paste0("slideseq_stackedbar_total100_clusternum",clusternum,"_",method,".pdf"),width=15,height=10)
print(p)
dev.off()
p = plot_celltype_barplot_each100(clusternum, celltypes, metadata_RCTD,method,cbp2,textsize=30)
pdf(paste0("slideseq_stackedbar_each100_clusternum",clusternum,"_",method,".pdf"),width=15,height=10)
print(p)
dev.off()

```

### Section 8: GSEA analysis.

Differential gene in regions.
```R
library(MAST)
library(fdrtool)
library(qvalue)

load("sim_SpatialPCA.RData") # load SpatialPCA slingshot object 
load(paste0("Puck_180430_6_zero_removed_nomt.rds"))
load("expr_0.05.RData")

sp_count = as.matrix(sp_count)
expr_counts = sp_count[,match(colnames(expr),colnames(sp_count))]
iCounts = expr_counts
lib_size <- colSums(iCounts)
tpm <- 1e+06*sweep(iCounts, 2, lib_size, "/")
tpm <- log2(tpm+1)
tpm = as.matrix(tpm)
diff_gene = list()
diff_gene_full = list()
for(celltypenum in 1:max(meta_data$SpatialPCA_Louvain)){

print(celltypenum)
cellType = meta_data$SpatialPCA_Louvain
cellType[which(cellType!=paste0(celltypenum))]="0"
cellType <- as.factor(cellType)
sca <- FromMatrix(tpm,  cData=data.frame(cellType=cellType))
freq_expressed <- 0.1
thres <- thresholdSCRNACountMatrix(assay(sca), nbins = 5, min_per_bin = 5, conditions = cellType) 
assays(sca) <- list(thresh=thres$counts_threshold, tpm=assay(sca))
expressed_genes <- freq(sca) > freq_expressed
sca <- sca[expressed_genes,]
ngeneson <- apply(iCounts,2,function(x) mean(x>0))
CD <- colData(sca)
CD$ngeneson <- ngeneson
CD$cngeneson <- CD$ngeneson-mean(ngeneson)
colData(sca) <- CD
fit <- zlm(~ cellType + cngeneson, sca = sca)
summaryDt <- summary(fit, doLRT=paste0("cellType",celltypenum))
summaryDt <- summaryDt$datatable
res.mast <- merge(summaryDt[contrast==paste0("cellType",celltypenum) & component=='H',.(primerid, `Pr(>Chisq)`)], summaryDt[contrast==paste0("cellType",celltypenum) & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid')
res.mast = res.mast[order(res.mast[,2]),]
diff_gene[[celltypenum]] = res.mast$primerid[which(res.mast[,2]<=0.05)]
diff_gene_full[[celltypenum]] = res.mast
}

for(celltypenum in 1:max(meta_data$SpatialPCA_Louvain)){
  pval = unlist(diff_gene_full[[celltypenum]][,2])
  fdr.hat.comp<-fdrtool(pval, statistic="pvalue")
  diff_gene[[celltypenum]] = diff_gene_full[[celltypenum]]$primerid[which(fdr.hat.comp$qval<=0.05)]
}

```

Use gene sets from the Molecular Signatures Database (MSigDB) available from the Broad Institute

```R
library(gprofiler2)
library(ggplot2)
library(tidyverse)
load("expr_0.05.RData")
load(paste0("Puck_180430_6_zero_removed_nomt.rds"))
SCTcount = sp_count[,match(colnames(expr),colnames(sp_count))]

h: hallmark gene sets hallmark gene sets as Gene Symbols  h.all.v7.4.symbols.gmt
c2: curated gene sets KEGG gene sets as Gene Symbols  c2.cp.kegg.v7.4.symbols.gmt
c4: computational gene sets cancer modules as Gene Symbols  c4.cm.v7.4.symbols.gmt
c5: Ontology gene sets all GO gene sets as Gene Symbols  c5.go.v7.4.symbols.gmt
c5: Ontology gene sets Human Phenotype Ontology as Gene Symbols  c5.hpo.v7.4.symbols.gmt
c7: immunologic signature gene sets all immunologic signature gene sets as Gene Symbols c7.all.v7.4.symbols.gmt

#-------------------
# make gmt files
#-------------------

# install.packages("ActivePathways")
library(ActivePathways)

gmt <- read.GMT("h.c255all.v7.4.symbols.gmt")
gmt_filter = gmt[c(1,2)]
count = 0
for(item in 1:length(gmt)){
  print(item)

  itemid = gmt[[item]]$id 
  if(length(grep("GSE",itemid))==0){
    count = count + 1
    gmt_filter[[count]]=gmt[[item]] 
  }
}

write.GMT(gmt_filter, "h.c255all.v7.4.symbols.filter.gmt")
upload_GMT_file(gmtfile = "h.c255all.v7.4.symbols.filter.gmt")

> upload_GMT_file(gmtfile = "h.c255all.v7.4.symbols.filter.gmt")
Your custom annotations ID is gp__Cb5w_84zs_UEE
You can use this ID as an 'organism' name in all the related enrichment tests against this custom source.
Just use: gost(my_genes, organism = 'gp__Cb5w_84zs_UEE')
[1] "gp__Cb5w_84zs_UEE"

```

Regionally differential gene.
```R
load("slideseq_diff_gene.RData")
datt = list()
for(num in 1:8){
  print(num)
topGOnum = 10
bg = rownames(SCTcounts)
gostres <- gost(query = diff_gene[[num]], 
                organism = "gp__Cb5w_84zs_UEE", ordered_query = FALSE, 
                multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "custom", custom_bg = bg, 
                numeric_ns = "", 
                #sources = c("GO:MF","GO:BP", "GO:CC","KEGG"), 
                as_short_link = FALSE)
gostres$result = gostres$result[order(gostres$result$p_value),]
gostres$result$Source = unlist(lapply(gostres$result$term_id,strsplit_func))
datt[[num]] = gostres$result[1:10,]
datt[[num]]$log10p = -log10(datt[[num]]$p_value)
}

spatial_domain = c("Granule inner sublayer","Molecular layer","Granule outer sublayer",
  "Purkinje layer","Choroid plexus","Granule middle sublayer","White matter","Cerebellar nucleus")
p=list()
for(num in 1:8){
  print(num)
p[[num]]=ggplot(data=datt[[num]], aes(x=fct_reorder(tolower(term_id),log10p), y=log10p, fill=Source,label = ifelse(significant ==TRUE, "*",""))) +
  geom_bar(stat="identity", position=position_dodge(),color="black",width=0.8)+
  #scale_fill_continuous(low='#F0E442', high='red', guide=guide_colorbar(reverse=TRUE)) + 
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
  labs(title=paste0(spatial_domain[num]),x="Biological terms", y = "-log10(p value)")+
  coord_flip()+
  theme_classic()+
  #geom_text(vjust = 0.5) +
  geom_text(vjust = 1, nudge_y = 0.5)+
  #ylim(0,1)+
    theme(plot.title = element_text(size = 20),
              text = element_text(size = 20,color="black"),
              #axis.title = element_text(face="bold"),
              axis.text.x=element_text(size = 20,color="black") ,
              legend.position = "right")# +
}

pdf(paste0("GSEA_SpatialPCA_msigdb_slideseq.pdf"),width=15,height=5)
for(num in 1:8){
  print(p[[num]])
}
dev.off()
```
Pseudotime related genes.
```R 

load("slideseq_pseudotime_gene_SpatialPCA.RData")
SCTcount = sp_count[,match(colnames(expr),colnames(sp_count))]
genelist = list()
geneset1 = as.character(data1$gene_traj1)[which(data1$p_value_traj1<0.05/dim(SCTcount)[1])]
geneset2 = as.character(data2$gene_traj2)[which(data2$p_value_traj2<0.05/dim(SCTcount)[1])]
geneset3 = as.character(data3$gene_traj3)[which(data3$p_value_traj3<0.05/dim(SCTcount)[1])]
genelist[[1]] = geneset1
genelist[[2]] = geneset2
genelist[[3]] = geneset3
load("/net/mulan/disk2/shanglu/Projects/spatialPCA/manuscript_v4/slideseq/pseudotime_related_geneset_granule.RData")
cellid = which(meta_data$SpatialPCA_Louvain %in% c("1","6","3"))
SCTcount_granule = SCTcount[,cellid]
data1 = na.omit(data.frame(gene_traj1, p_value_traj1))
geneset_granule = as.character(data1$gene_traj1)[which(data1$p_value_traj1<0.05/dim(SCTcount_granule)[1])]
genelist[[4]] = geneset_granule

titles = c(paste0("Pseudotime trajectory ",1:3),"Pseudotime trajectory - Granule layer")

datt = list()
for(num in 1:4){
  print(num)
topGOnum = 10
bg = rownames(SCTcounts)
gostres <- gost(query = genelist[[num]], 
                organism = "gp__Cb5w_84zs_UEE", ordered_query = FALSE, 
                multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "custom", custom_bg = bg, 
                numeric_ns = "", 
                #sources = c("GO:MF","GO:BP", "GO:CC","KEGG"), 
                as_short_link = FALSE)
gostres$result = gostres$result[order(gostres$result$p_value),]
gostres$result$Source = unlist(lapply(gostres$result$term_id,strsplit_func))
datt[[num]] = gostres$result[1:10,]
datt[[num]]$log10p = -log10(datt[[num]]$p_value)
}

p=list()
for(num in 1:4){
  print(num)
p[[num]]=ggplot(data=datt[[num]], aes(x=fct_reorder(tolower(term_id),log10p), y=log10p, fill=Source,label = ifelse(significant ==TRUE, "*",""))) +
  geom_bar(stat="identity", position=position_dodge(),color="black",width=0.8)+
  #scale_fill_continuous(low='#F0E442', high='red', guide=guide_colorbar(reverse=TRUE)) + 
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
  labs(title=paste0(titles[num]),x="Biological terms", y = "-log10(p value)")+
  coord_flip()+
  theme_classic()+
  #geom_text(vjust = 0.5) +
  geom_text(vjust = 1, nudge_y = 0.5)+
  #ylim(0,1)+
    theme(plot.title = element_text(size = 20),
              text = element_text(size = 20,color="black"),
              #axis.title = element_text(face="bold"),
              axis.text.x=element_text(size = 20,color="black") ,
              legend.position = "right")# +
}

pdf(paste0("GSEA_SpatialPCA_pseudotime_msigdb_slideseq.pdf"),width=15,height=5)
for(num in 1:4){
  print(p[[num]])
}
dev.off()


```








