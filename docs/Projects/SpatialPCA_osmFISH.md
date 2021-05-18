---
layout: default
title: osmFISH analysis
nav_order: 3
has_children: false
parent: SpatialPCA
permalink: /docs/Projects/SpatialPCA/osmFISH
---


##### We illustrate the benefits of SpatialPCA through four different downstream analyses: spatial transcriptomics visualization, spatial domain detection, spatial trajectory inference on the tissue, and high-resolution spatial map reconstruction. 


##### osmFISH expressin and coordinate data are downloaded from [https://github.com/linnarsson-lab/osmFISH_celltype_analysis/blob/master/osmFISH_SScortex_mouse_all_cells.loom](https://github.com/linnarsson-lab/osmFISH_celltype_analysis/blob/master/osmFISH_SScortex_mouse_all_cells.loom) and processed into txt files.

### Section 1: Data processing.

```R
  library(Seurat)
  library(data.table)
  library(ggplot2)
 osm_exprs = read.table(file = paste0('osmFISH_expression.txt'))
 osm_locs = read.table(file = paste0('osmFISH_cell_coordinates.txt'))
 osm_locs = osm_locs[rownames(osm_locs) %in% colnames(osm_exprs),]
 osmfish = CreateSeuratObject(counts = osm_exprs,meta.data=osm_locs,min.cells = 10, min.features = 10)
 osmfish = SCTransform(osmfish, 
               return.only.var.genes = FALSE, 
               variable.features.n = NULL, 
               variable.features.rv.th = 1.3)
 SCTcounts = osmfish@assays$SCT@counts
 expr = osmfish@assays$SCT@scale.data
 info = scale(osm_locs[match(colnames(expr), rownames(osm_locs)),])
 save(expr, info, file = "osmFISH_expr.RData")

```


### Section 2: Run SpatialPCA.

SpatialPCA:
```R
load("osmFISH_expr.RData")
ls()
## [1] "expr" "info"
expr=scale_expr(expr)
dat = data_prepare_func(expr, info)
bandwidth = bandwidth_select(expr, info,method="Silverman")
K=kernel_build(kernelpara="gaussian", dat$ED2,bandwidth) 

# Set the number of PCs. Here we pre-calculate some matrices that will be useful in SpatialPCA functions for large data. In this way we could save some time in each iteration.
PC_num = 20
dat$YMt = t(dat$YM)
dat$KYM = K%*%dat$YMt
dat$K = K
eigen_res = eigs_sym(K, k=PC_num, which = "LM")
dat$delta = eigen_res$values
dat$U = eigen_res$vectors

Est_para = SpatialPCA_estimate_parameter_largedata(dat_input=dat,PCnum=20)
Est_W = SpatialPCA_estimate_W_largedata(Est_para$par, dat,PCnum=20)
Est_Z = SpatialPCA_estimate_Z_largedata(Est_para$par,dat,Est_W,PCnum=20)


```

PCA and NMF:
```R
dat$Z_pca=get_PCA(expr, PCnum=20)
dat$Z_NMF = get_NMF(expr, PCnum=20)
```

HMRF:
```python
# prepare files in r
setwd("~/osmFISH/HMRF")
ids = c(1:dim(expr)[2])
exprdata = as.data.frame(t(expr))
writeexpr = cbind(ids,exprdata)
write.table(writeexpr, file = "~/osmFISH/HMRF/expr_marker.txt", row.names=F,col.names=F, quote=F )
writelocation = cbind(ids,0,info )
write.table(writelocation, file = "~/osmFISH/HMRF/location.txt", row.names=F,col.names=F, quote=F )
gene = rownames(expr)
write.table(gene, file = "~/osmFISH/HMRF/genes", row.names=F,col.names=F, quote=F )

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

directory = "~/osmFISH/HMRF"
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
os.chdir('~/osmFISH/HMRF') 
this_hmrf = HMRFInstance("osmFISH", outdir, new_dset, 12,  (0, 0.5, 30), tolerance=1e-10)
this_hmrf.init(nstart=1000, seed=-1)
this_hmrf.run()
stop = timeit.default_timer()
print('Time: ', stop - start)  
visualize.domain(this_hmrf, 12, 10, dot_size=8, size_factor=10, outfile="cluster12.visualize.beta.%.1f.png" % 5.0)


# collect results
k = 12
name = "osmFISH"
betas = seq(from=0,to=14.5,by=0.5)
python_path = "/net/mulan/home/shanglu/anaconda3/bin/python3.7"
get_result_path = system.file("python", "get_result2.py", package = 'Giotto')
output_data = "~/osmFISH/HMRF/spatial.HMRF.cluster"
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
p_PCs = plot_factor_value(info,dat$Z_spatial,textmethod="SpatialPCA",pointsize=1.5,textsize=15)
ggarrange(p_PCs[[1]], p_PCs[[2]], p_PCs[[3]],
          ncol = 3, nrow = 1)

p_PCs = plot_factor_value(info,dat$Z_pca,textmethod="PCA",pointsize=1.5,textsize=15)
ggarrange(p_PCs[[1]], p_PCs[[2]], p_PCs[[3]],
          ncol = 3, nrow = 1)

p_PCs = plot_factor_value(info,dat$Z_NMF,textmethod="NMF",pointsize=1.5,textsize=15)
ggarrange(p_PCs[[1]], p_PCs[[2]], p_PCs[[3]],
          ncol = 3, nrow = 1)

```
Visualization by RGB plots.
```R
library(ggpubr)

p1 = plot_RGB_tSNE(info,dat$Z_spatial)
p2 = plot_RGB_tSNE(info,dat$Z_pca)
p3 = plot_RGB_tSNE(info,dat$Z_NMF)
pdf("RGB_tSNE.pdf",width=15, height=5)
ggarrange(p1, p2, p3, 
           labels = c("SpatialPCA", "PCA", "NMF"),
          ncol = 3, nrow = 1)
dev.off()

p1 = plot_RGB_UMAP(info,dat$Z_spatial)
p2 = plot_RGB_UMAP(info,dat$Z_pca)
p3 = plot_RGB_UMAP(info,dat$Z_NMF)
pdf("RGB_UMAP.pdf",width=15, height=5)
ggarrange(p1, p2, p3, 
          labels = c("SpatialPCA", "PCA", "NMF"),
          ncol = 3, nrow = 1)
dev.off()

```

### Section 4: Use clustering results from SpatialPCA to obtain spatial domains.
We used walktrap clustering results in the main analysis. Louvain clustering results are in the supplementary figures. We tried a sequence of number of nearest neighbors in the graph construction for both Walktrap method and the Louvain method. We selected the number of nearest neighbors that could result in 12 clusters based on original osmFISH paper. When multiple number of nearest neighbors could result in same number of clusters, we usually choose the largest number among them.
```R
walktrap_cluster_SpatialPCA = walktrap_clustering(knearest = seq(from=20,to=150,by=1), latent_dat=dat$Z_spatial)
walktrap_cluster_PCA = walktrap_clustering(knearest = seq(from=20,to=150,by=1), latent_dat=dat$Z_pca)
walktrap_cluster_NMF = walktrap_clustering(knearest = seq(from=20,to=150,by=1), latent_dat=dat$Z_NMF)


meta_data = as.data.frame(info)
meta_data$HMRF = cluster_HMRF[[21]]
meta_data$SpatialPCA_Walktrap = walktrap_cluster_SpatialPCA$cluster_label[[63]]
meta_data$PCA_Walktrap = walktrap_cluster_PCA$cluster_label[[14]]
meta_data$NMF_Walktrap = walktrap_cluster_NMF$cluster_label[[1]]

```

Compare with ground truth annotated by pathologist through ARI, NMI, and Pseudo-R2.
```R

# use original annotation from osmFISH paper
df_coordinates = read.csv("~/osmFISH/df_coordinates", sep="\t")
df_coordinates = t(df_coordinates)
colnames(df_coordinates) = df_coordinates[1,]
df_coordinates1 = df_coordinates[-1,] 
df_coordinates1 = as.data.frame(df_coordinates1)
cellid = paste0("cell_",df_coordinates1$CellID)
df_coordinates1$cellid = cellid
df_coordinates_use = df_coordinates1[match(colnames(expr), cellid),]
tmp = info[match(colnames(expr),rownames(info)),]
meta = cbind(tmp, info)
meta = cbind(meta,df_coordinates_use )


meta_data = cbind(meta_data, meta)
save(meta_data, file = "osmfish_meta_data.RData")

library(DescTools)
library(nnet)
library(mcca)

labels = metadata$Region
NMFfactor = as.data.frame(t(dat$Z_NMF))
colnames(NMFfactor) = paste0("NMFfactor",1:20)
spatialPCAfactor = as.data.frame(t(dat$Z_spatial))
colnames(spatialPCAfactor) = paste0("spatialPCAfactor",1:20)
PCAfactor = as.data.frame(t(dat$Z_pca))
colnames(PCAfactor) = paste0("spatialPCAfactor",1:20)

ind=which(as.character(labels)=="Excluded")
labels_update = labels[-ind]    

fit <- multinom(labels_update ~ as.matrix(spatialPCAfactor)[-ind ,], maxit = 1000, MaxNWts = 2000,model = TRUE)
SpatialPCA_R2 = PseudoR2(fit,c("McFaddenAdj"))

fit <- multinom(labels_update ~ as.matrix(PCAfactor)[-ind ,], maxit = 1000, MaxNWts = 2000,model = TRUE)
PCA_R2 = PseudoR2(fit,c("McFaddenAdj"))

fit <- multinom(labels_update ~ as.matrix(NMFfactor)[-ind ,], maxit = 1000, MaxNWts = 2000,model = TRUE)
NMF_R2 = PseudoR2(fit,c("McFaddenAdj"))

fit <- multinom(labels_update ~ as.factor(meta_data$HMRF[-ind]), maxit = 1000, MaxNWts = 2000,model = TRUE)
HMRF_label_R2 = PseudoR2(fit,c("McFaddenAdj"))


SpatialPCA_R2
PCA_R2
NMF_R2
HMRF_label_R2

#> SpatialPCA_R2
#McFaddenAdj 
#  0.9505127 
#> PCA_R2
#McFaddenAdj 
#  0.4780278 
#> NMF_R2
#McFaddenAdj 
#   0.532203 
#> HMRF_label_R2
#McFaddenAdj 
#  0.6119262 

R_square = c(SpatialPCA_R2, PCA_R2, NMF_R2, HMRF_label_R2)
method = c("SpatialPCA", "PCA","NMF","HMRF")

datt = data.frame(R_square, method)
datt$method = factor(datt$method, levels=method, order=T)
pdf(paste0("osmFISH_PseudoR2.pdf"),width=5,height=5)
ggplot(data=datt, aes(x=method, y=R_square, fill=method)) +
  geom_bar(stat="identity", position=position_dodge(),color="black",width=0.8)+
  scale_fill_brewer(palette="Paired")+
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
  labs(title="Pseudo R2",x="", y = "")+
  theme(legend.position="bottom") +
  theme_classic()+ylim(0,1)+
    theme(plot.title = element_text(size = 22),
              text = element_text(size = 22),
              #axis.title = element_text(face="bold"),
              #axis.text.x=element_text(size = 20) ,
              legend.position = "right")# +
dev.off()


ind=which(as.character(labels)=="Excluded")
compare_table = data.frame(
              meta_data$SpatialPCA_Walktrap, 
              meta_data$PCA_Walktrap,
              meta_data$NMF_Walktrap, 
              meta_data$HMRF
              )

compare_table$cluster_true = labels
compare_table = compare_table[-ind,]
my_nmi = c()
my_acc = c()
for(count in 1:4){
  cluster_result = as.factor(compare_table[,count])
  cluster_true = as.factor(compare_table[,5])
my_nmi[count] = compare(cluster_result,cluster_true, method = "nmi" )
my_acc[count] = adjustedRandIndex(cluster_result, cluster_true)

}
my_nmi
my_acc

#> my_nmi
#[1] 0.7504445 0.3402852 0.4075793 0.6281221 
#> my_acc
#[1] 0.6874798 0.2191418 0.3064690 0.5729922 
NMI = my_nmi
ARI = my_acc
method = c("SpatialPCA", "PCA","NMF","HMRF")
library(wesanderson)
datt = data.frame(NMI, ARI,method)
datt$method = factor(datt$method, levels=method, order=T)
pdf(paste0("osmfish_NMI.pdf"),width=5,height=5)
ggplot(data=datt, aes(x=method, y=NMI, fill=method)) +
  geom_bar(stat="identity", position=position_dodge(),color="black",width=0.8)+
  scale_fill_brewer(palette="Paired")+
  labs(title="NMI",x="", y = "")+
  theme(legend.position="bottom") +
  theme_classic()+
  ylim(0,1)+
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
  theme(plot.title = element_text(size = 22),
              text = element_text(size = 22),
              #axis.title = element_text(face="bold"),
              #axis.text.x=element_text(size = 20) ,
              legend.position = "right")# +
dev.off()


datt = data.frame(NMI, ARI,method)
datt$method = factor(datt$method, levels=method, order=T)
pdf(paste0("osmfish_ARI.pdf"),width=5,height=5)
ggplot(data=datt, aes(x=method, y=ARI, fill=method)) +
  geom_bar(stat="identity", position=position_dodge(),color="black",width=0.8)+
  scale_fill_brewer(palette="Paired")+
  labs(title="ARI",x="", y = "")+
  theme(legend.position="bottom") +
  theme_classic()+
  ylim(0,1)+
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
  theme(plot.title = element_text(size = 22),
              text = element_text(size = 22),
              #axis.title = element_text(face="bold"),
              #axis.text.x=element_text(size = 20) ,
              legend.position = "right")# +
dev.off()



```

### Section 5: Clustering results visualization.

```R
loc1 = info[,1]
loc2 = info[,2]

# SpatialPCA walktrap
cbp_SpatialPCA_Walktrap <- c( "chartreuse3", "coral2", "tan1",  "yellow",  "dodgerblue",
  "#CC79A7",  "steelblue",   "cadetblue1","skyblue2","#52854C","#FFDB6D","slateblue")
pdf("osmfish_Plot_12_clusters_SpatialPCA_Walktrap.pdf",width=5,height=10)
plot_cluster(loc1,loc2,meta_data$SpatialPCA_Walktrap,pointsize=3,text_size=22 ,"",cbp_SpatialPCA_Walktrap)
dev.off()


# PCA walktrap
cbp_PCA_Walktrap <- c( "#CC79A7", "coral2", "dodgerblue",  "yellow",  "tan1",
  "#52854C",  "steelblue",   "cadetblue1","#FFDB6D","chartreuse3","skyblue2","slateblue")
pdf("osmfish_Plot_12_clusters_PCA_Walktrap.pdf",width=5,height=10)
plot_cluster(loc1,loc2,meta_data$PCA_Walktrap,pointsize=3,text_size=22 ,"",cbp_PCA_Walktrap)
dev.off()

# NMF walktrap
cbp_NMF_Walktrap <- c( "slateblue", "#CC79A7", "skyblue2",  "tan1",  "yellow",
  "coral2",  "chartreuse3",   "#52854C","cadetblue1","#FFDB6D","steelblue","dodgerblue")
pdf("osmfish_Plot_12_clusters_NMF_Walktrap.pdf",width=5,height=10)
plot_cluster(loc1,loc2,meta_data$NMF_Walktrap,pointsize=3,text_size=22 ,"",cbp_NMF_Walktrap)
dev.off()


# HMRF
cbp_HMRF <- c( "chartreuse3", "cadetblue1", "#FFDB6D",  "slateblue",  "yellow",
  "tan1",  "skyblue",   "coral2","steelblue","dodgerblue","#CC79A7","skyblue2")
pdf("osmfish_Plot_12_clusters_HMRF.pdf",width=5,height=10)
plot_cluster(loc1,loc2,meta_data$HMRF,pointsize=3,text_size=22 ,"",cbp_HMRF)
dev.off()


```


### Section 6: High resolution map reconstruction.
```R
Z_high = high_resolution(info, K, kernelpara="gaussian",ED=dat$ED, est_log_tau = dat$Est_para$par,est_W = Est_W[[1]] ,est_sigma0 = Est_W[[2]][1,1],est_Z = Est_Z,PCnum=20)
walktrap_SpatialPCA_highresolution = walktrap_clustering(knearest = 560, Z_high$Z_star)
cbp_highresolution <- c( "#52854C", "tan1", "steelblue",  "yellow",  "coral2",
  "dodgerblue",  "#FFDB6D",   "skyblue2","cadetblue1","slateblue",
  "chartreuse3","#CC79A7")
loc1=unlist(Z_high$Location_star[,1])
loc2=unlist(Z_high$Location_star[,2])
Z_high_clusters = walktrap_SpatialPCA_highresolution$cluster_label[[1]]
p3 = plot_cluster(loc1,loc2,Z_high_clusters,pointsize=0.5,text_size=10 ,title_in="High-resolution",cbp_highresolution)
print(p3)
```

### Section 7: Trajectory analysis.
In original data:
```R
library(slingshot)
library(SingleCellExperiment)

method = "SpatialPCA"
sim<- SingleCellExperiment(assays = list(counts = as.matrix(dat$Z_spatial)))
reducedDims(sim) <- SimpleList(DRM = t(dat$Z_spatial))
colData(sim)$Walktrap <- metadata$SpatialPCA_Walktrap_para1    
sim  <-slingshot(sim, clusterLabels = 'Walktrap', reducedDim = 'DRM',start.clus="12" )
summary(sim@colData@listData)
#> summary(sim@colData@listData)
#                  Length Class  Mode     
#Walktrap          5275   -none- character
#slingPseudotime_1 5275   -none- numeric  
#slingPseudotime_2 5275   -none- numeric  
#slingPseudotime_3 5275   -none- numeric  

pseudotime_traj1 = sim@colData@listData$slingPseudotime_1
pseudotime_traj2 = sim@colData@listData$slingPseudotime_2
pseudotime_traj3 = sim@colData@listData$slingPseudotime_3
clusterlabels = meta_data$SpatialPCA_Walktrap
gridnum = 10
color_in = c(  "plum1", "dodgerblue","mediumaquamarine",  "palegreen4", "chocolate1","lightblue2","#F0E442",  "black","#CC79A7","mediumpurple","seagreen1")
p_traj1 = plot_trajectory(pseudotime_traj1, info,clusterlabels,gridnum,color_in,pointsize=2 ,arrowlength=0.2,arrowsize=0.8,textsize=15 )
p_traj2 = plot_trajectory(pseudotime_traj2, info,clusterlabels,gridnum,color_in,pointsize=2 ,arrowlength=0.2,arrowsize=0.8,textsize=15 )
p_traj3 = plot_trajectory(pseudotime_traj3, info,clusterlabels,gridnum,color_in,pointsize=2 ,arrowlength=0.2,arrowsize=0.8,textsize=15 )
print(ggarrange( p_traj1[[4]],p_traj2[[4]],p_traj3[[4]],p_traj1[[1]],p_traj2[[1]],p_traj3[[1]],
          ncol = 3, nrow = 2))
```
High-resolution spatial map reconstruction:
```R
simhigh <- SingleCellExperiment(assays = list(counts = as.matrix(Z_high$Z_star)))
reducedDims(simhigh) <- SimpleList(DRM = t(Z_high$Z_star))
colData(simhigh)$louvain <- factor(louvain_SpatialPCA_highresolution_cluster)    
simhigh  <-slingshot(simhigh, clusterLabels = 'louvain', reducedDim = 'DRM',start.clus="7" )
summary(simhigh@colData@listData)

#> summary(simhigh@colData@listData)
#                  Length Class  Mode   
#louvain           21100  factor numeric
#slingPseudotime_1 21100  -none- numeric
#slingPseudotime_2 21100  -none- numeric
#slingPseudotime_3 21100  -none- numeric
#slingPseudotime_4 21100  -none- numeric


Z_high_pseudotime_traj1 = simhigh@colData@listData$slingPseudotime_1
Z_high_pseudotime_traj2 = simhigh@colData@listData$slingPseudotime_2
Z_high_pseudotime_traj3 = simhigh@colData@listData$slingPseudotime_3
Z_high_pseudotime_traj4 = simhigh@colData@listData$slingPseudotime_4

cluster = Z_high_clusters
gridnum = 20
color_in = c(  "plum1", "palegreen4","mediumaquamarine",  "dodgerblue", "chocolate1",
            "#F0E442","lightblue2",  "black")
p_traj1 = plot_trajectory(Z_high_pseudotime_traj1, Z_high$Location_star,cluster,gridnum,color_in,pointsize=0.5 ,arrowlength=0.2,arrowsize=0.8,textsize=15 )
p_traj2 = plot_trajectory(Z_high_pseudotime_traj2, Z_high$Location_star,cluster,gridnum,color_in,pointsize=0.5 ,arrowlength=0.2,arrowsize=0.8,textsize=15 )
p_traj3 = plot_trajectory(Z_high_pseudotime_traj3, Z_high$Location_star,cluster,gridnum,color_in,pointsize=0.5 ,arrowlength=0.2,arrowsize=0.8,textsize=15 )
p_traj4 = plot_trajectory(Z_high_pseudotime_traj4, Z_high$Location_star,cluster,gridnum,color_in,pointsize=0.5 ,arrowlength=0.2,arrowsize=0.8,textsize=15 )
print(ggarrange( p_traj1[[4]],p_traj2[[4]],p_traj3[[4]],p_traj4[[4]],p_traj1[[1]],p_traj2[[1]],p_traj3[[1]],p_traj4[[1]],
          ncol = 4, nrow = 2))
```

### Section 8: Region annotations from osmFISH paper distribution in each spatial domain from SpatialPCA.
```R
ind=which(as.character(labels)=="Excluded")
meta_data_filter = meta_data[-ind,]

cbp1 <- c( "steelblue", "slateblue",  "dodgerblue","cadetblue1",  "#CC79A7",
 "#FFDB6D",  "coral2",  "skyblue2", "tan1","yellow",
  "chartreuse3")
clusternum = 12
regionnum = 11
meta_data_filter$Region = as.factor(meta_data_filter$Region)
percentage = matrix(0,clusternum,regionnum)
for(k in 1:clusternum){
metadata_sub = meta_data_filter[which(meta_data_filter$SpatialPCA_Walktrap==paste0(k) ),]
match_type = metadata_sub$Region
percentage[k,] = round(unlist(table(match_type))/dim(metadata_sub)[1]*100,2)
}


rownames(percentage) = paste0("Cluster",1:clusternum)
Region = names(table(match_type))
colnames(percentage) = Region
percentage_vec = c(percentage)
cluster_vec = c(rep(c(paste0("Cluster",1:clusternum)),length(regionnum)))
Region = c(rep(Region,each=clusternum))
datt = data.frame(cluster_vec, percentage_vec,Region)
datt$cluster_vec = factor(cluster_vec, level=paste0("Cluster",1:regionnum))
datt$Region = as.character(datt$Region)
pdf("osmfish_region_each100_Stack_barplot_all_SpatialPCA_Walktrap.pdf",width=10,height=4)
  ggplot(datt, aes(y = percentage_vec,
             x = factor(cluster_vec ), fill = Region)) +        ## global aes
  scale_fill_manual(values=cbp1)+
  geom_bar(position="stack", stat="identity",width=0.7,color="black") +
  theme_void()+xlab("")+ylab("")+ggtitle("SpatialPCA") +
  theme(plot.title = element_text(size = 20),
              text = element_text(size = 20),
              #axis.title = element_text(face="bold"),
              #axis.text.x=element_text(size = 12,angle = 60,hjust = 1) ,
              #axis.text.x=element_blank(),
              legend.position = "right")# +
dev.off()


clusternum = 12
regionnum = 11
meta_data_filter$Region = as.factor(meta_data_filter$Region)
percentage = matrix(0,clusternum,regionnum)
for(k in 1:clusternum){
metadata_sub = meta_data_filter[which(meta_data_filter$SpatialPCA_Walktrap==paste0(k) ),]
match_type = metadata_sub$Region
percentage[k,] = round(unlist(table(match_type))/dim(meta_data_filter)[1]*100,2)
}

rownames(percentage) = paste0("Cluster",1:clusternum)
Region = names(table(match_type))
colnames(percentage) = Region
percentage_vec = c(percentage)
cluster_vec = c(rep(c(paste0("Cluster",1:clusternum)),length(regionnum)))
Region = c(rep(Region,each=clusternum))
datt = data.frame(cluster_vec, percentage_vec,Region)
datt$cluster_vec = factor(cluster_vec, level=paste0("Cluster",1:clusternum))
datt$Region = as.character(datt$Region)
pdf("osmfish_region_total100_Stack_barplot_all_SpatialPCA_Walktrap.pdf",width=10,height=4)
  ggplot(datt, aes(y = percentage_vec,
             x = factor(cluster_vec ), fill = Region)) +        ## global aes
  scale_fill_manual(values=cbp1)+
  geom_bar(position="stack", stat="identity",width=0.7,color="black") +
  theme_classic()+xlab("")+ylab("")+
  theme(plot.title = element_text(size = 20),
              text = element_text(size = 20),
              #axis.title = element_text(face="bold"),
              #axis.text.x=element_text(size = 12,angle = 60,hjust = 1) ,
              #axis.text.x=element_blank(),
              legend.position = "right")# +
dev.off()

```

### Section 9: Cell type annotations from osmFISH paper distribution in each spatial domain from SpatialPCA.
```R

cbp1 <- c( "chartreuse3", "tan1", "yellow",  "#FFDB6D", "steelblue", 
   "skyblue2", "#CC79A7",   "cadetblue1","coral2","slateblue",
  "dodgerblue","#BAD968", "#E1989D", "#D48F90", "#F1BA63", "#C8C4D3", "#E1F3B8", "#F78377",
"#FCCDE5","#D6CA65", "#E2B379" ,"#C6D989" ,"#CAADC4" ,"#86B1CD" ,"#E0DEC5",
"#8DD3C7" ,"#B4B2A3","#A6A1B4" ,"#B7E3BF" ,"#E1D3B7", "#F8F7B7")
# col=colorRampPalette(brewer.pal(8, "Set3"))(20)[sample(20,replace=F)]

set.seed=1234
cbp_stack <- c(colorRampPalette(brewer.pal(8, "Set3"))(11),colorRampPalette(brewer.pal(8, "Set2"))(11),colorRampPalette(brewer.pal(8, "Set1"))(11))[sample(31,replace=F)]

# SpatialPCA
clusternum = 12
ClusterNamenum = 31
meta_data_filter$ClusterName = as.factor(meta_data_filter$ClusterName)
percentage = matrix(0,clusternum,ClusterNamenum)
for(k in 1:clusternum){
metadata_sub = meta_data_filter[which(meta_data_filter$SpatialPCA_Walktrap==k ),]
match_type = metadata_sub$ClusterName
percentage[k,] = round(unlist(table(match_type))/dim(metadata_sub)[1]*100,2)
}

rownames(percentage) = paste0("Cluster",1:clusternum)
ClusterName = as.character(unique(meta_data_filter$ClusterName))
colnames(percentage) = ClusterName
rownames(percentage) = paste0("Cluster",1:clusternum)
percentage_vec = c(percentage)
cluster_vec = c(rep(c(paste0("Cluster",1:clusternum)),length(ClusterName)))
ClusterName = c(rep(ClusterName,each=clusternum))
datt = data.frame(cluster_vec, percentage_vec,ClusterName)
datt$cluster_vec = factor(cluster_vec, level=paste0("Cluster",1:clusternum))
datt$ClusterName = as.character(datt$ClusterName)
pdf("osmfish_ClusterName_each100_Stack_barplot_all_SpatialPCA_Walktrap.pdf",width=10,height=6)
  ggplot(datt, aes(y = percentage_vec,
             x = cluster_vec, fill = ClusterName)) +        ## global aes
  scale_fill_manual(values=cbp_stack)+
  geom_bar(position="stack", stat="identity",width=0.7,color="black") +
  theme_void()+xlab("")+ylab("")+
  theme(plot.title = element_text(size = 15),
              text = element_text(size = 15),
              #axis.title = element_text(face="bold"),
              #axis.text.x=element_text(size = 12,angle = 60,hjust = 1) ,
              #axis.text.x=element_blank(),
              legend.position = "right")# +
dev.off()


clusternum = 12
ClusterNamenum = 31
meta_data_filter$ClusterName = as.factor(meta_data_filter$ClusterName)
percentage = matrix(0,clusternum,ClusterNamenum)
for(k in 1:clusternum){
metadata_sub = meta_data_filter[which(meta_data_filter$SpatialPCA_Walktrap==k ),]
match_type = metadata_sub$ClusterName
percentage[k,] = round(unlist(table(match_type))/dim(meta_data_filter)[1]*100,2)
}

rownames(percentage) = paste0("Cluster",1:clusternum)
ClusterName = as.character(unique(meta_data_filter$ClusterName))
colnames(percentage) = ClusterName
rownames(percentage) = paste0("Cluster",1:clusternum)
percentage_vec = c(percentage)
cluster_vec = c(rep(c(paste0("Cluster",1:clusternum)),length(ClusterName)))
ClusterName = c(rep(ClusterName,each=clusternum))
datt = data.frame(cluster_vec, percentage_vec,ClusterName)
datt$cluster_vec = factor(cluster_vec, level=paste0("Cluster",1:clusternum))
datt$ClusterName = as.character(datt$ClusterName)
pdf("osmfish_ClusterName_total100_Stack_barplot_all_SpatialPCA_Walktrap.pdf",width=10,height=6)
  ggplot(datt, aes(y = percentage_vec,
             x = cluster_vec, fill = ClusterName)) +        ## global aes
  scale_fill_manual(values=cbp_stack)+
  geom_bar(position="stack", stat="identity",width=0.7,color="black") +
  theme_classic()+xlab("")+ylab("")+
  theme(plot.title = element_text(size = 15),
              text = element_text(size = 15),
              #axis.title = element_text(face="bold"),
              #axis.text.x=element_text(size = 12,angle = 60,hjust = 1) ,
              #axis.text.x=element_blank(),
              legend.position = "right")# +
dev.off()


```




