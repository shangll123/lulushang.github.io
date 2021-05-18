---
layout: default
title: ST data analysis
nav_order: 2
has_children: false
parent: SpatialPCA
permalink: /docs/Projects/SpatialPCA/ST
---


##### We illustrate the benefits of SpatialPCA through four different downstream analyses: spatial transcriptomics visualization, spatial domain detection, spatial trajectory inference on the tissue, and high-resolution spatial map reconstruction. 

##### ST breast tumor data are downloaded from [https://github.com/almaan/her2st](https://github.com/almaan/her2st).

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

Following data processing codes from [https://github.com/almaan/her2st](https://github.com/almaan/her2st):

```R
meta_data <- read.csv("../res/ST-cluster/sample.csv",header=TRUE, stringsAsFactors=FALSE,sep=",")
meta_data$patient_id = c()
for(i in 1:dim(meta_data)[1]){
  meta_data$patient_id[i] = paste0(meta_data$patient[i],meta_data$cluster[i] )
}
rownames(meta_data) = meta_data$patient_id
samples <- list.files(pattern = ".tsv", path = "../data/ST-cnts/", full.names = T)
names(samples) <- substr(do.call(rbind, strsplit(samples, split = "/"))[, 5], start = 1, stop = 2)
imgs <- list.files(path = "../data/ST-imgs/", recursive = T, full.names = T, pattern = ".jpg")
names(imgs) <- do.call(rbind, strsplit(imgs, split = "/"))[, 6]
ids <- names(samples)
infoTable <- data.frame(samples, imgs = imgs[ids], ids, patient_id = substr(x = ids, start = 1, stop = 1), stringsAsFactors = FALSE)

tmp = meta_data[match(infoTable$ids, meta_data$patient_id),]

infoTable <- cbind(infoTable, tmp)
infoTable[, -c(11:28)]

#Subset infoTable to include specified datasets.
infoTable$spotfiles <- list.files(path = "../data/ST-spotfiles", full.names = T)[1:36]
head(infoTable)

## Load data
#Load all patient datasets and merge into one Seurat object per patient. Each gene has to bre present in at least 20 spots per sample and each spot has to have at least 300 unique features (genes).
seu.list <- lapply(unique(infoTable$patient_id), function(s) {
    InputFromTable(infotable = subset(infoTable, patient_id == s), 
                      min.gene.spots = 20,
                      min.spot.feature.count = 300,
                      platform = "1k")
}
)

# remove ring genes
seu.list <- lapply(seu.list, function(seu) {
  subset(seu, features = setdiff(rownames(seu@assays$RNA@counts), ring.genes))
})


#Calculate some QC metrics
total.qc <- do.call(rbind, lapply(seu.list, function(se) {
  data.frame(total_UMIs = sum(se@assays$RNA@counts), nSpots = ncol(se))
}))


# QC

qcMat <- do.call(rbind, lapply(1:length(seu.list), function(i) {
    seu <- seu.list[[i]]
    do.call(rbind, lapply(unique(seu[["ids", drop = T]]), function(id) {
        repMat <- seu@assays$RNA@counts[, seu[["ids", drop = T]] == id]
        nUMI <- Matrix::colSums(repMat)
        nGene <- apply(repMat, 2, function(x) sum(x > 0))
        data.frame(sample = id, 
                   avg.nUMI = round(mean(nUMI)),
                   median.nUMI = median(nUMI),
                   max.nUMI = max(nUMI),
                   min.nUMI = min(nUMI),
                   avg.nGene = round(mean(nGene)),
                   median.nGene = median(nGene),
                   min.nGene = min(nGene),
                   max.nGene = max(nGene),
                   nSpots = ncol(repMat))
    }))
}))

qcMat

```

Prepare count matrix and normalized matrix. We used H1 sample in the manuscript.
```R


# default variable.features.rv.th is 1.3, original https://github.com/almaan/her2st used 1.1
seu.list.1.3 = seu.list
seu.list.1.3 <- lapply(seu.list.1.3, function(seu) {
  SCTransform(seu, 
              vars.to.regress = c("ids"), 
              return.only.var.genes = FALSE, 
              variable.features.n = NULL, 
              variable.features.rv.th = 1.3)
})

for(num in 1:length(seu.list.1.3)){
  print(num)
  seu.list.single = seu.list.1.3[[num]]
  save(seu.list.single, file = paste0("~/her2st/data/seu.list.single",num,".rv1.3.RData"))
}

her2stdatanum = 0
for(num in 1:8){
  load(paste0("~/her2st/data/seu.list.single",num,".rv1.3.RData"))
  count = 0
  for(id in unique(seu.list.single[["ids", drop = T]])){
    count = count + 1
    her2stdatanum = her2stdatanum + 1
    print(her2stdatanum)
    SCTcounts = seu.list.single@assays$SCT@counts[, seu.list.single[["ids", drop = T]] == id]
    SCTscaled = seu.list.single@assays$SCT@scale.data[, seu.list.single[["ids", drop = T]] == id]
    ind = which(seu.list.single@tools$Staffli@meta.data$sample == count)
    metadata = seu.list.single@tools$Staffli@meta.data[ind,]
    print(dim(metadata))
    save(SCTcounts, SCTscaled, metadata, file = paste0("~/her2st/data/her2stdata",her2stdatanum,".rv1.3.RData"))
  }
}
```

We use SPARK to select spatial genes as input for SpatialPCA.
```R

her2stdatanum=34 # this is the H1 sample
load(paste0("~/her2st/data/her2stdata",her2stdatanum,".rv1.3.RData"))

num=8 # H1 sample belongs to the 8th individual.
load(paste0("~/her2st/data/seu.list.single",num,".rv1.3.RData"))

H_count = seu.list.single@assays$RNA@counts
rawcount = H_count[match(rownames(SCTcounts),rownames(H_count)),match(colnames(SCTcounts),colnames(H_count))]
# > dim(rawcount)
# [1] 10053   607

library(SPARK)
info = scale(metadata[,c(3:4)])
info = as.data.frame(info)
info[,2]=-info[,2] # this is to flip the figure by y axis for easier visualization. 
rownames(info) = colnames(rawcount)
spark <- CreateSPARKObject(counts=rawcount, 
                                 location=info[,1:2],
                                 percentage = 0.1, 
                                 min_total_counts = 10)

spark@lib_size <- apply(rawcount, 2, sum)
spark <- spark.vc(spark, 
          covariates = NULL, 
                    lib_size = spark@lib_size, 
                    num_core = 5,
                    verbose = F,
                    fit.model="gaussian")
spark <- spark.test(spark, 
                    check_positive = T, 
                    verbose = T)

ress = spark@res_mtest[,c("combined_pvalue","adjusted_pvalue")]  

sum(ress$adjusted_pvalue<0.1)
sum(ress$adjusted_pvalue<0.05)
sum(ress$adjusted_pvalue<0.01)
# > sum(ress$adjusted_pvalue<0.1)
# [1] 368
# > sum(ress$adjusted_pvalue<0.05)
# [1] 319
# > sum(ress$adjusted_pvalue<0.01)
# [1] 237

# we explored different p value cut offs for spatial gene selection in the sensitivity analysis.
expr = SCTscaled[match(rownames(ress)[which(ress$adjusted_pvalue<=0.05)], rownames(SCTscaled)),]
save(expr, info,file = "expr_0.05.RData")
expr = SCTscaled[match(rownames(ress)[which(ress$adjusted_pvalue<=0.01)], rownames(SCTscaled)),]
save(expr, info, file = "expr_0.01.RData")
expr = SCTscaled[match(rownames(ress)[which(ress$adjusted_pvalue<=0.1)], rownames(SCTscaled)),]
save(expr,info,  file = "expr_0.1.RData")

```

### Section 2: Run SpatialPCA.

SpatialPCA:
```R
load("expr_0.05.RData")
ls()
## [1] "expr" "info"
expr=scale_expr(expr)
dat = data_prepare_func(expr, info)
bandwidth = bandwidth_select(expr, info,method="SJ")
K=kernel_build(kernelpara="gaussian", dat$ED2,bandwidth) 
Est_para = SpatialPCA_estimate_parameter(dat_input=dat,PCnum=20)
Est_W = SpatialPCA_estimate_W(Est_para$par, dat,PCnum=20)
Est_Z = SpatialPCA_estimate_Z(Est_para$par,dat,Est_W,PCnum=20)

```

PCA and NMF:
```R
dat$Z_pca=get_PCA(expr, PCnum=20)
dat$Z_NMF = get_NMF(expr, PCnum=20)
```

HMRF:
```python
# prepare files in r
setwd("~/her2st/HMRF")
ids = c(1:dim(expr)[2])
exprdata = as.data.frame(t(expr))
writeexpr = cbind(ids,exprdata)
write.table(writeexpr, file = "~/her2st/HMRF/expr_marker.txt", row.names=F,col.names=F, quote=F )
writelocation = cbind(ids,0,info )
write.table(writelocation, file = "~/her2st/HMRF/location.txt", row.names=F,col.names=F, quote=F )
gene = rownames(expr)
write.table(gene, file = "~/her2st/HMRF/genes", row.names=F,col.names=F, quote=F )

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

directory = "~/her2st/HMRF"
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
os.chdir('~/her2st/HMRF') 
this_hmrf = HMRFInstance("ST", outdir, new_dset, 7,  (0, 0.5, 30), tolerance=1e-10)
this_hmrf.init(nstart=1000, seed=-1)
this_hmrf.run()
stop = timeit.default_timer()
print('Time: ', stop - start)  
visualize.domain(this_hmrf, 7, 10, dot_size=8, size_factor=10, outfile="cluster7.visualize.beta.%.1f.png" % 5.0)


# collect results
k = 7
name = "ST"
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
p_PCs = plot_factor_value(info,dat$Z_spatial,textmethod="SpatialPCA",pointsize=3.5,textsize=15)
ggarrange(p_PCs[[1]], p_PCs[[2]], p_PCs[[3]],
          ncol = 3, nrow = 1)

p_PCs = plot_factor_value(info,dat$Z_pca,textmethod="PCA",pointsize=3.5,textsize=15)
ggarrange(p_PCs[[1]], p_PCs[[2]], p_PCs[[3]],
          ncol = 3, nrow = 1)

p_PCs = plot_factor_value(info,dat$Z_NMF,textmethod="NMF",pointsize=3.5,textsize=15)
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
We used walktrap clustering results in the main analysis. Louvain clustering results are in the supplementary figures. We tried a sequence of number of nearest neighbors in the graph construction for both Walktrap method and the Louvain method. We selected the number of nearest neighbors that could result in 7 clusters based on histological image obtained by pathologists. When multiple number of nearest neighbors could result in same number of clusters, we usually choose the largest number among them.
```R
walktrap_cluster_SpatialPCA = walktrap_clustering(knearest = c(seq(from=10,to=50,by=1)), latent_dat=dat$Z_spatial)
walktrap_cluster_PCA = walktrap_clustering(knearest = c(seq(from=10,to=50,by=1)), latent_dat=dat$Z_pca)
walktrap_cluster_NMF = walktrap_clustering(knearest = c(seq(from=10,to=50,by=1)), latent_dat=dat$Z_NMF)

louvain_cluster_SpatialPCA = louvain_clustering(knearest = c(seq(from=10,to=50,by=1)), latent_dat=dat$Z_spatial)
louvain_cluster_PCA = louvain_clustering(knearest = c(seq(from=10,to=50,by=1)), latent_dat=dat$Z_pca)
louvain_cluster_NMF = louvain_clustering(knearest = c(seq(from=10,to=50,by=1)), latent_dat=dat$Z_NMF)


her2stdatanum=34
load(paste0("~/her2st/data/her2stdata",her2stdatanum,".rv1.3.RData"))
path_anno = "~/her2st/her2st-master/data/ST-pat/lbl"
metadata$order = c(1:dim(metadata)[1])
anno = read.csv(paste0(path_anno,"/H1_labeled_coordinates.tsv"),sep="\t")
colnames(anno) = c("Row.names"  , "adj_x" ,"adj_y",  "pixel_x"  ,"pixel_y" ,"label")
myanno = merge(metadata, anno, by=c("adj_x","adj_y"))
anno_label = myanno[order(myanno$order),]
labels = anno_label$label

meta_data = cbind(info,anno_label)
meta_data$SpatialPCA_Walktrap = walktrap_cluster_SpatialPCA$cluster_label[[22]]
meta_data$SpatialPCA_Louvain = louvain_cluster_SpatialPCA$cluster_label[[12]]
meta_data$PCA_Walktrap = walktrap_cluster_PCA$cluster_label[[25]]
meta_data$PCA_Louvain = louvain_cluster_PCA$cluster_label[[10]]
meta_data$NMF_Walktrap = walktrap_cluster_NMF$cluster_label[[32]]
meta_data$NMF_Louvain = louvain_cluster_NMF$cluster_label[[30]]
meta_data$HMRF = cluster_HMRF[[21]]
```

Compare with ground truth annotated by pathologist through ARI, NMI, and Pseudo-R2.
```R

library(DescTools)
#https://www.rdocumentation.org/packages/DescTools/versions/0.99.37/topics/PseudoR2
library(nnet)
library(mcca)

labels = anno_label$label
NMFfactor = as.data.frame(t(dat$Z_NMF))
colnames(NMFfactor) = paste0("NMFfactor",1:20)
spatialPCAfactor = as.data.frame(t(dat$Z_spatial))
colnames(spatialPCAfactor) = paste0("spatialPCAfactor",1:20)
PCAfactor = as.data.frame(t(dat$Z_pca))
colnames(PCAfactor) = paste0("spatialPCAfactor",1:20)

ind=which(as.character(labels)=="undetermined")
labels_update = labels[-ind]    

fit <- multinom(labels_update ~ as.matrix(spatialPCAfactor)[-ind ,], maxit = 1000, MaxNWts = 2000,model = TRUE)
SpatialPCA_R2 = PseudoR2(fit,c("McFaddenAdj"))

fit <- multinom(labels_update ~ as.matrix(PCAfactor)[-ind ,], maxit = 1000, MaxNWts = 2000,model = TRUE)
PCA_R2 = PseudoR2(fit,c("McFaddenAdj"))

fit <- multinom(labels_update ~ as.matrix(NMFfactor)[-ind ,], maxit = 1000, MaxNWts = 2000,model = TRUE)
NMF_R2 = PseudoR2(fit,c("McFaddenAdj"))

fit <- multinom(labels_update ~ as.factor(meta_data$HMRF[-ind]), maxit = 1000, MaxNWts = 2000,model = TRUE)
HMRF_label_R2 = PseudoR2(fit,c("McFaddenAdj"))

NMF_R2
PCA_R2
SpatialPCA_R2
HMRF_label_R2

# > NMF_R2
# McFaddenAdj 
#   0.5889822 
# > PCA_R2
# McFaddenAdj 
#   0.6029148 
# > SpatialPCA_R2
# McFaddenAdj 
#   0.7763146 
# > HMRF_label_R2
# McFaddenAdj 
#   0.4669912 

R_square = c(SpatialPCA_R2, PCA_R2, NMF_R2, HMRF_label_R2)
method = c("SpatialPCA", "PCA","NMF","HMRF")

datt = data.frame(R_square, method)
datt$method = factor(datt$method, levels=method, order=T)
pdf(paste0("PseudoR2.pdf"),width=5,height=5)
ggplot(data=datt, aes(x=method, y=R_square, fill=method)) +
  geom_bar(stat="identity", position=position_dodge(),color="black",width=0.8)+
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
  labs(title="PseudoR2",x="", y = "")+
  theme(legend.position="bottom") +
  theme_classic()+ylim(0,1)+
    theme(plot.title = element_text(size = 15),
              text = element_text(size = 15),
              #axis.title = element_text(face="bold"),
              #axis.text.x=element_text(size = 20) ,
              legend.position = "right")# +
dev.off()


ind=which(as.character(labels)=="undetermined")
compare_table = data.frame(meta_data$HMRF, 
              meta_data$SpatialPCA_Walktrap, #meta_data$SpatialPCA_Louvain, 
              meta_data$PCA_Walktrap,#meta_data$PCA_Louvain,
              meta_data$NMF_Walktrap #meta_data$NMF_Louvain
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

> my_nmi
[1] 0.4904179 0.5461827 0.4056973 0.3472605
> my_acc
[1] 0.3874706 0.4327208 0.3301797 0.2561279


NMI = c(0.5461827,0.4056973 ,0.3472605,0.4904179)
ARI = c(0.4327208,0.3301797 ,0.2561279,0.3874706)

method = c("SpatialPCA", "PCA","NMF","HMRF")

library(wesanderson)
datt = data.frame(NMI, ARI,method)
datt$method = factor(datt$method, levels=method, order=T)
pdf(paste0("NMI.pdf"),width=5,height=5)
ggplot(data=datt, aes(x=method, y=NMI, fill=method)) +
  geom_bar(stat="identity", position=position_dodge(),color="black",width=0.8)+
  scale_fill_brewer(palette="Paired")+
  labs(title="NMI",x="", y = "")+
  theme(legend.position="bottom") +
  theme_classic()+
  ylim(0,1)+
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
  theme(plot.title = element_text(size = 15),
              text = element_text(size = 15),
              #axis.title = element_text(face="bold"),
              #axis.text.x=element_text(size = 20) ,
              legend.position = "right")# +
dev.off()


datt = data.frame(NMI, ARI,method)
datt$method = factor(datt$method, levels=method, order=T)
pdf(paste0("ARI.pdf"),width=5,height=5)
ggplot(data=datt, aes(x=method, y=ARI, fill=method)) +
  geom_bar(stat="identity", position=position_dodge(),color="black",width=0.8)+
  scale_fill_brewer(palette="Paired")+
  labs(title="ARI",x="", y = "")+
  theme(legend.position="bottom") +
  theme_classic()+
  ylim(0,1)+
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
  theme(plot.title = element_text(size = 15),
              text = element_text(size = 15),
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
cbp1 <- c(  "plum1", "dodgerblue","mediumaquamarine",  "palegreen4", "chocolate1",
            "lightblue2","#F0E442",  "red","#CC79A7","mediumpurple","seagreen1")
pdf("Plot_7_clusters_SpatialPCA_Walktrap.pdf",width=5,height=5)
plot_cluster(loc1,loc2,meta_data$SpatialPCA_Walktrap,pointsize=4.5,text_size=10 ,"",cbp1)
dev.off()

# PCA walktrap
cbp1 <- c(  "plum1", "chocolate1","dodgerblue",  "palegreen4",  "lightblue2",
           "#F0E442","mediumaquamarine",  "red","#CC79A7","mediumpurple","seagreen1")
pdf("Plot_7_clusters_PCA_Walktrap.pdf",width=5,height=5)
plot_cluster(loc1,loc2,meta_data$PCA_Walktrap,pointsize=4.5,text_size=10 ,"",cbp1)
dev.off()

# NMF walktrap
cbp1 <- c(  "plum1", "dodgerblue",  "#F0E442","chocolate1", "mediumaquamarine", 
            "lightblue2","palegreen4", "red","#CC79A7","mediumpurple","seagreen1")
pdf("Plot_7_clusters_NMF_Walktrap.pdf",width=5,height=5)
plot_cluster(loc1,loc2,meta_data$NMF_Walktrap,pointsize=4.5,text_size=10 ,"",cbp1)
dev.off()

# HMRF
cbp1 <- c(  "chocolate1", "dodgerblue","#F0E442",  "palegreen4", "mediumaquamarine",
            "plum1","lightblue2",  "red","#CC79A7","mediumpurple","seagreen1")
pdf("Plot_7_clusters_HMRF.pdf",width=5,height=5)
plot_cluster(loc1,loc2,meta_data$HMRF,pointsize=4.5,text_size=10 ,"",cbp1)
dev.off()

```


### Section 6: High resolution map reconstruction.
```R
Z_high = high_resolution(info, K, kernelpara="gaussian",ED=dat$ED, est_log_tau = dat$Est_para$par,est_W = Est_W[[1]] ,est_sigma0 = Est_W[[2]][1,1],est_Z = Est_Z,PCnum=20)
walktrap_SpatialPCA_highresolution = walktrap_clustering(knearest = 84, Z_high$Z_star)
cbp_high = c(  "plum1", "palegreen4","mediumaquamarine",  "dodgerblue", "chocolate1", "#F0E442","lightblue2")
loc1=unlist(Z_high$Location_star[,1])
loc2=unlist(Z_high$Location_star[,2])
Z_high_clusters = walktrap_SpatialPCA_highresolution$cluster_label[[1]]
p3 = plot_cluster(loc1,loc2,Z_high_clusters,pointsize=2,text_size=10 ,title_in="High-resolution",cbp_high)
print(p3)
```

### Section 7: Trajectory analysis.
In original data:
```R
sim<- SingleCellExperiment(assays = expr)
reducedDims(sim) <- SimpleList(DRM = t(dat$Z_spatial))
colData(sim)$Walktrap <- factor(meta_data$SpatialPCA_Walktrap)    
sim  <-slingshot(sim, clusterLabels = 'Walktrap', reducedDim = 'DRM',start.clus="5" ) # define starting point in cancer region
summary(sim@colData@listData)

pseudotime_traj1 = sim@colData@listData$slingPseudotime_1
pseudotime_traj2 = sim@colData@listData$slingPseudotime_2
pseudotime_traj3 = sim@colData@listData$slingPseudotime_3
clusterlabels = meta_data$SpatialPCA_Walktrap
gridnum = 10
color_in = c(  "plum1", "dodgerblue","mediumaquamarine",  "palegreen4", "chocolate1","lightblue2","#F0E442",  "black","#CC79A7","mediumpurple","seagreen1")
p_traj1 = plot_trajectory(pseudotime_traj1, info,clusterlabels,gridnum,color_in,pointsize=4 ,arrowlength=0.2,arrowsize=0.8,textsize=15 )
p_traj2 = plot_trajectory(pseudotime_traj2, info,clusterlabels,gridnum,color_in,pointsize=4 ,arrowlength=0.2,arrowsize=0.8,textsize=15 )
p_traj3 = plot_trajectory(pseudotime_traj3, info,clusterlabels,gridnum,color_in,pointsize=4 ,arrowlength=0.2,arrowsize=0.8,textsize=15 )
print(ggarrange( p_traj1[[4]],p_traj2[[4]],p_traj3[[4]],p_traj1[[1]],p_traj2[[1]],p_traj3[[1]],
          ncol = 3, nrow = 2))
```
High-resolution spatial map reconstruction:
```R
simhigh <- SingleCellExperiment(Z_high$Z_star)
reducedDims(simhigh) <- SimpleList(DRM = t(Z_high$Z_star))
colData(simhigh)$Walktrap <- factor(Z_high_clusters)    
simhigh  <-slingshot(sim, clusterLabels = 'Walktrap', reducedDim = 'DRM',start.clus="5" ) 
summary(simhigh@colData@listData)

Z_high_pseudotime_traj1 = simhigh@colData@listData$slingPseudotime_1
Z_high_pseudotime_traj2 = simhigh@colData@listData$slingPseudotime_2
Z_high_pseudotime_traj3 = simhigh@colData@listData$slingPseudotime_3
cluster = Z_high_clusters
gridnum = 20
color_in = c(  "plum1", "palegreen4","mediumaquamarine",  "dodgerblue", "chocolate1",
            "#F0E442","lightblue2",  "black")
p_traj1 = plot_trajectory(Z_high_pseudotime_traj1, Z_high$Location_star,cluster,gridnum,color_in,pointsize=1.5 ,arrowlength=0.2,arrowsize=0.8,textsize=15 )
p_traj2 = plot_trajectory(Z_high_pseudotime_traj2, Z_high$Location_star,cluster,gridnum,color_in,pointsize=1.5 ,arrowlength=0.2,arrowsize=0.8,textsize=15 )
p_traj3 = plot_trajectory(Z_high_pseudotime_traj3, Z_high$Location_star,cluster,gridnum,color_in,pointsize=1.5 ,arrowlength=0.2,arrowsize=0.8,textsize=15 )

print(ggarrange( p_traj1[[4]],p_traj2[[4]],p_traj3[[4]],p_traj1[[1]],p_traj2[[1]],p_traj3[[1]],
          ncol = 3, nrow = 2))
```

### Section 8: Cell type deconvolution.
We used RCTD to perform cell type deconvolution. The reference data is from [Single-Cell Map of Diverse Immune Phenotypesin the Breast Tumor Microenvironment](https://www.sciencedirect.com/science/article/pii/S0092867418307232?via%3Dihub)
```R
# devtools::install_github("dmcable/RCTD", build_vignettes = TRUE)
library(RCTD)
library(Matrix)
# GSE114725
# - Azizi, Elham, Ambrose J. Carr, George Plitas, Andrew E. Cornish, Catherine Konopacki, Sandhya Prabhakaran, Juozas Nainys, et al. “[Single-Cell Map of Diverse Immune Phenotypes in the Breast Tumor Microenvironment](https://doi.org/10.1016/j.cell.2018.05.060).” Cell, June 2018. - scRNA-seq of immune cells in BRCA - continuous activation of T cells, no macrophage polarization. inDrop and 10X platforms. 47,016 CD45+ cells from 8 primary breast carcinomas. 83 clusters, tested by cross-validation. [Data 1](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114727), [Data 2](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114725), [Data 3](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114724)
GSE114725 = read.csv("~/her2st/her2st_reference/raw_corrected.csv",header=T)
# they already removed clusters that were not used in their analysis
countmat = GSE114725[,-c(1:5)]
rownames(countmat) = paste0(GSE114725$patient, GSE114725$tissue,GSE114725$cellid)
meta_data = GSE114725[,1:5]
meta_data$sample = rownames(countmat)
countmat=as.matrix(countmat)
countmat = t(countmat)

raw.data = countmat
meta_data$nUMI = colSums(raw.data)
rownames(meta_data) = as.character(meta_data$sample)
meta_data$liger_ident_coarse = as.factor(meta_data$cluster)
reference = Seurat::CreateSeuratObject(raw.data, meta.data = meta_data)
saveRDS(reference, paste(getwd(),"SCRef_her2st_dana.RDS",sep="/"))
dref <- create_downsampled_data(reference, getwd())
dref <- create_downsampled_data(dref, getwd(), n_samples = 300)
create_downsampled_data(dref, getwd(), n_samples = 25)

library(doParallel)
myRCTD <- create.RCTD(puck, reference, max_cores = 5)
myRCTD <- run.RCTD(myRCTD, doublet_mode = TRUE)
results <- myRCTD@results
norm_weights = sweep(results$weights, 1, rowSums(results$weights), '/') 
cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
spatialRNA <- myRCTD@spatialRNA
doublets <- results$results_df[results$results_df$spot_class == "doublet_certain",] 
doub_occur <- table(doublets$second_type, doublets$first_type) 


cluster_names = c( "1" ,  "2" ,  "3" ,  "4" ,  "5" ,  "6"  , "8"  , "9" ,
 "11" , "12" , "15" ,"16" , "19" , "20" ,
  "21" , "23"  ,"24" , "25" , "27", "28" ,
  "29"  ,"30" , "31" , "35" , "36" , "37" , "38" , "39" , "40",  
  "41" , "42" , "43" , "44" , "45" , "46" , "47" , "48" , "49" , "50" , 
  "51" ,"53",  "54" , "55"  ,"56" , "57" , "58" , "59" , "60" , 
  "61" , "62" , "63" , "64" , "65" , "66" , "67" , "68" , "70" , 
  "71" , "72" , "73","74" , "75" , "76",  "77" , "78" , "80" ,
   "81"  ,"82" , "83" , "84" , "85" , "86" , "87" , "88" , "89" , "90" , 
   "91" , "92" , "93" , "94", "95"  )

cluster_annotation = c("T_CD4","T_CD4","T_CD8","NK","T_CD8","B","B","NKT",
  "T_CD8","NK","T_CD8","T_CD4","T_CD4","T_CD4","mDC","MACROPHAGE","MONOCYTE_precursor","MACROPHAGE","MONOCYTE_precursor","MACROPHAGE",
  "NEUTROPHIL","NK","MAST","NK","MONOCYTE","mDC","mDC","B","MONOCYTE",
  "pDC","MONOCYTE_precursor","T_CD8","MAST","T_CD8","T_Reg","T_CD4","B","MONOCYTE_precursor","NK",
  "T_CD8","B","MONOCYTE","B","T_Reg","B","T_CD8","T_CD8","T_CD4",
  "B","T_CD8","NK","MONOCYTE_precursor","B","T_CD4","MONOCYTE","MONOCYTE","T_CD4",
  "T_CD4","T_CD8","MAST","T_CD8","T_CD4","T_CD8","T_CD8","NK","T_Reg",
  "T_CD4","T_Reg","T_CD8","MONOCYTE","T_CD8","MONOCYTE","T_Reg","NK","NEUTROPHIL","B",
  "MONOCYTE","NK","T_CD4","MONOCYTE","T_CD4")

anno = data.frame(cluster_names,cluster_annotation)
celltypes = as.character(results$results_df$first_type)
result_table = results$results_df
result_table$clusternum_type1 = as.character(result_table$first_type)
result_table$match_type1 = anno$cluster_annotation[match(result_table$clusternum_type1, anno$cluster_names)]
result_table$clusternum_type2 = as.character(result_table$second_type)
result_table$match_type2 = anno$cluster_annotation[match(result_table$clusternum_type2, anno$cluster_names)]

metadata = cbind(metadata,result_table)
meta_data$match_type = metadata$match_type

```
Plot cell types in each cluster.
```R
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(ggthemes)

cbp2 = c("#FFDB6D", "#C4961A", "#F4EDCA", "tomato","#C3D7A4",  "#4E84C4", "#293352","#52854C",
"#D16103", "deepskyblue1", "cadetblue3","lightblue1","plum1","chartreuse3")

percentage = matrix(0,7,13)
for(k in 1:7){
metadata_sub = meta_data[which(meta_data$SpatialPCA_Walktrap==k ),]
match_type = metadata_sub$match_type
percentage[k,] = round(unlist(table(match_type))/dim(metadata_sub)[1]*100,2)
}
rownames(percentage) = paste0("Cluster",1:7)
celltype = c("B cell","Macrophage","Mast cell", "mDC","Monocyte","Monocyte precursor","Neutrophil","NK cell","NKT cell","pDC cell","T CD4", "T CD8","T Reg")
colnames(percentage) = celltype
rownames(percentage) = paste0("Cluster",1:7)
percentage_vec = c(percentage)
cluster_vec = c(rep(c(paste0("Cluster",1:7)),length(celltype)))
CellType = c(rep(celltype,each=7))
datt = data.frame(cluster_vec, percentage_vec,CellType)
datt$cluster_vec = factor(cluster_vec, level=paste0("Cluster",1:7))

pdf("Stack_barplot_all.pdf",width=8,height=8)
  ggplot(datt, aes(y = percentage_vec,
             x = factor(cluster_vec ), fill = CellType)) +        ## global aes
  scale_fill_manual(values=cbp2)+
  geom_bar(position="stack", stat="identity",width=0.7) +
  theme_void()+xlab("")+ylab("")+
  theme(plot.title = element_text(size = 15),
              text = element_text(size = 15),
              #axis.title = element_text(face="bold"),
              #axis.text.x=element_text(size = 12,angle = 60,hjust = 1) ,
              #axis.text.x=element_blank(),
              legend.position = "bottom")# +
dev.off()


percentage = matrix(0,7,13)
for(k in 1:7){
metadata_sub = meta_data[which(meta_data$SpatialPCA_Walktrap==k ),]
match_type = metadata_sub$match_type
percentage[k,] = round(unlist(table(match_type))/dim(metadata)[1]*100,2)
}
rownames(percentage) = paste0("Cluster",1:7)
celltype = c("B cell","Macrophage","Mast cell", "mDC","Monocyte","Monocyte precursor","Neutrophil","NK cell","NKT cell","pDC cell","T CD4", "T CD8","T Reg")
colnames(percentage) = celltype
rownames(percentage) = paste0("Cluster",1:7)
percentage_vec = c(percentage)
cluster_vec = c(rep(c(paste0("Cluster",1:7)),length(celltype)))
CellType = c(rep(celltype,each=7))
datt = data.frame(cluster_vec, percentage_vec,CellType)
datt$cluster_vec = factor(cluster_vec, level=paste0("Cluster",1:7))

pdf("Stack_barplot_each.pdf",width=8,height=8)
  ggplot(datt, aes(y = percentage_vec,
             x = factor(cluster_vec ), fill = CellType)) +        ## global aes
  scale_fill_manual(values=cbp2)+
  geom_bar(position="stack", stat="identity",color="black",width=0.7) +
  theme_classic()+xlab("")+ylab("")+
  theme(plot.title = element_text(size = 15),
              text = element_text(size = 15),
              #axis.title = element_text(face="bold"),
              #axis.text.x=element_text(size = 12,angle = 60,hjust = 1) ,
              #axis.text.x=element_blank(),
              legend.position = "bottom")# +
dev.off()

```

### Section 9: GSEA analysis.

Differential gene in regions.
```R

library(MAST)
library(fdrtool)
library(qvalue)


her2stdatanum=34
load(paste0("~/her2st/data/her2stdata",her2stdatanum,".rv1.3.RData"))
num=8
load(paste0("~/her2st/data/seu.list.single",num,".rv1.3.RData"))
H_count = seu.list.single@assays$RNA@counts
rawcount = H_count[match(rownames(SCTcounts),rownames(H_count)),match(colnames(SCTcounts),colnames(H_count))]

iCounts = rawcount
lib_size <- colSums(iCounts)
tpm <- 1e+06*sweep(iCounts, 2, lib_size, "/")
tpm <- log2(tpm+1)
tpm = as.matrix(tpm)


diff_gene_full = list()
for(celltypenum in 1:max(meta_data$SpatialPCA_Walktrap)){

print(celltypenum)
cellType = meta_data$SpatialPCA_Walktrap
cellType[which(cellType!=paste0(celltypenum))]="0"
cellType <- as.factor(cellType)
sca <- FromMatrix(tpm,  cData=data.frame(cellType=cellType))
# sca <- FromMatrix(tpm,  cData=data.frame(cellType=cellType, batch=batch))
freq_expressed <- 0.1
thres <- thresholdSCRNACountMatrix(assay(sca), nbins = 10, min_per_bin = 50, conditions = cellType) # worked, current results, also default
assays(sca) <- list(thresh=thres$counts_threshold, tpm=assay(sca))
expressed_genes <- freq(sca) > freq_expressed
sca <- sca[expressed_genes,]
ngeneson <- apply(iCounts,2,function(x) mean(x>0))
CD <- colData(sca)
CD$ngeneson <- ngeneson
CD$cngeneson <- CD$ngeneson-mean(ngeneson)
colData(sca) <- CD
## differential expression
fit <- zlm(~ cellType + cngeneson, sca = sca)
summaryDt <- summary(fit, doLRT=paste0("cellType",celltypenum))
summaryDt <- summaryDt$datatable
res.mast <- merge(summaryDt[contrast==paste0("cellType",celltypenum) & component=='H',.(primerid, `Pr(>Chisq)`)], summaryDt[contrast==paste0("cellType",celltypenum) & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid')
res.mast = res.mast[order(res.mast[,2]),]
diff_gene_full[[celltypenum]] = res.mast
}


for(celltypenum in 1:max(meta_data$SpatialPCA_Walktrap)){
  pval = unlist(diff_gene_full[[celltypenum]][,2])
  fdr.hat.comp<-fdrtool(pval, statistic="pvalue")
  diff_gene[[celltypenum]] = diff_gene_full[[celltypenum]]$primerid[which(fdr.hat.comp$qval<=0.05)]
}

```

Pseudotime related genes.
```R
p_value_traj1 = c()
gene_traj1 = c()
p_value_traj2 = c()
gene_traj2 = c()
p_value_traj3 = c()
gene_traj3 = c()
for(i in 1:dim(rawcount)[1]){
  print(i)
    expr = rawcount[i,]
    dat = data.frame(expr, sim@colData@listData$slingPseudotime_1)
    colnames(dat) = c("expr","pseudotime")
    dat = na.omit(dat)
    res = try(summary(m1 <- glm.nb(expr ~ pseudotime, data = dat)))
    if(isTRUE(class(res)=="try-error")) { next } else { 
      p_value_traj1[i] =  res$coefficients[2,4] 
      gene_traj1[i] = rownames(rawcount)[i]
    } 

    dat = data.frame(expr, sim@colData@listData$slingPseudotime_2)
    colnames(dat) = c("expr","pseudotime")
    dat = na.omit(dat)
    res = try(summary(m1 <- glm.nb(expr ~ pseudotime, data = dat)))
    if(isTRUE(class(res)=="try-error")) { next } else { 
      p_value_traj2[i] =  res$coefficients[2,4] 
      gene_traj2[i] = rownames(rawcount)[i]
    } 

    dat = data.frame(expr, sim@colData@listData$slingPseudotime_3)
    colnames(dat) = c("expr","pseudotime")
    dat = na.omit(dat)
    res = try(summary(m1 <- glm.nb(expr ~ pseudotime, data = dat)))
    if(isTRUE(class(res)=="try-error")) { next } else { 
      p_value_traj3[i] =  res$coefficients[2,4] 
      gene_traj3[i] = rownames(rawcount)[i]
    } 
}


save(gene_traj1,gene_traj2,gene_traj3,p_value_traj1,p_value_traj2,p_value_traj3, file = "pseudotime_related_geneset.RData")

data1 = na.omit(data.frame(gene_traj1, p_value_traj1))
data2 = na.omit(data.frame(gene_traj2, p_value_traj2))
data3 = na.omit(data.frame(gene_traj3, p_value_traj3))

geneset1 = as.character(data1$gene_traj1)[which(data1$p_value_traj1<0.05/10053)]
geneset2 = as.character(data2$gene_traj2)[which(data2$p_value_traj2<0.05/10053)]
geneset3 = as.character(data3$gene_traj3)[which(data3$p_value_traj3<0.05/10053)]

write.table(geneset1, file = "pseudotime_geneset1.txt",col.names=F, row.names=F, quote=F)
write.table(geneset2, file = "pseudotime_geneset2.txt",col.names=F, row.names=F, quote=F)
write.table(geneset3, file = "pseudotime_geneset3.txt",col.names=F, row.names=F, quote=F)
write.table(rownames(rawcount), file = "pseudotime_allgene.txt",col.names=F, row.names=F, quote=F)


```
Use gene sets from the Molecular Signatures Database (MSigDB) available from the Broad Institute

```R

library(gprofiler2)
library(ggplot2)
library(tidyverse)

h: hallmark gene sets hallmark gene sets as Gene Symbols  h.all.v7.4.symbols.gmt
c2: curated gene sets KEGG gene sets as Gene Symbols  c2.cp.kegg.v7.4.symbols.gmt
c4: computational gene sets cancer modules as Gene Symbols  c4.cm.v7.4.symbols.gmt
c5: Ontology gene sets all GO gene sets as Gene Symbols  c5.go.v7.4.symbols.gmt
c5: Ontology gene sets Human Phenotype Ontology as Gene Symbols  c5.hpo.v7.4.symbols.gmt
c7: immunologic signature gene sets all immunologic signature gene sets as Gene Symbols c7.all.v7.4.symbols.gmt

# make gmt files
# install.packages("ActivePathways")
library(ActivePathways)
gmt <- read.GMT("h.c24557all.v7.4.symbols.gmt") # made this file by combining gmt files

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
write.GMT(gmt_filter, "h.c24557all.v7.4.symbols.filter.gmt")
upload_GMT_file(gmtfile = "h.c24557all.v7.4.symbols.filter.gmt")
> upload_GMT_file(gmtfile = "h.c24557all.v7.4.symbols.filter.gmt")
Your custom annotations ID is gp__iBGE_LwEo_qO4
You can use this ID as an 'organism' name in all the related enrichment tests against this custom source.
Just use: gost(my_genes, organism = 'gp__iBGE_LwEo_qO4')
[1] "gp__iBGE_LwEo_qO4"


# differential gene
strsplit_func = function(input){
  strsplit(input, split = "_")[[1]][1]
}

datt = list()
for(num in 1:7){
  print(num)
topGOnum = 10
bg = rownames(SCTcounts)
gostres <- gost(query = diff_gene[[num]], 
                organism = "gp__iBGE_LwEo_qO4", ordered_query = FALSE, 
                multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "custom", custom_bg = bg, 
                numeric_ns = "", 
                #sources = c("GO:MF","GO:BP", "GO:CC","KEGG"), 
                as_short_link = FALSE)
gostres$result = gostres$result[order(gostres$result$p_value),]
gostres$result$Source = unlist(lapply(gostres$result$term_id,strsplit_func))
gostres$result$Source[which(gostres$result$Source=="MODULE")]="CANCER MODULE"
datt[[num]] = gostres$result[1:10,]
datt[[num]]$log10p = -log10(datt[[num]]$p_value)
}

spatial_domain = c("Tumor surrounding region","Fibrous tissue near tumor","Fibrous tissue near normal glands",
  "Normal glands","Tumor region","Fat tissue","Immune region")
p=list()
for(num in 1:7){
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
pdf(paste0("GSEA_SpatialPCA_msigdb_her2st.pdf"),width=15,height=5)
for(num in 1:7){
  print(p[[num]])
}
dev.off()


# pseudotime related gene
datt = list()
for(num in 1:3){
  print(num)
topGOnum = 10
genelist = as.character(read.table(paste0("~/her2st/pseudotime_geneset",num,".txt"),header=F)$V1)
bg = rownames(SCTcounts)
gostres <- gost(query = genelist, 
                organism = "gp__iBGE_LwEo_qO4", ordered_query = FALSE, 
                multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "custom", custom_bg = bg, 
                numeric_ns = "", 
                #sources = c("GO:MF","GO:BP", "GO:CC","KEGG"), 
                as_short_link = FALSE)
gostres$result = gostres$result[order(gostres$result$p_value),]
gostres$result$Source = unlist(lapply(gostres$result$term_id,strsplit_func))
gostres$result$Source[which(gostres$result$Source=="MODULE")]="CANCER MODULE"
datt[[num]] = gostres$result[1:10,]
datt[[num]]$log10p = -log10(datt[[num]]$p_value)
}

p=list()
for(num in 1:3){
  print(num)
p[[num]]=ggplot(data=datt[[num]], aes(x=fct_reorder(tolower(term_id),log10p), y=log10p, fill=Source,label = ifelse(significant ==TRUE, "*",""))) +
  geom_bar(stat="identity", position=position_dodge(),color="black",width=0.8)+
  #scale_fill_continuous(low='#F0E442', high='red', guide=guide_colorbar(reverse=TRUE)) + 
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
  labs(title=paste0("Pseudotime trajectory ",num),x="Biological terms", y = "-log10(p value)")+
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

pdf(paste0("GSEA_SpatialPCA_pseudotime_msigdb_her2st.pdf"),width=12,height=5)
for(num in 1:3){
  print(p[[num]])
}
dev.off()

```

### Section 10: TLS region identification.
TLS gene expression heatmap
```R

mapDrugToColor<-function(annotations){
    colorsVector = ifelse(annotations["category"]=="Cluster1", 
        "red", ifelse(annotations["category"]=="Cluster2", 
        "orange",ifelse(annotations["category"]=="Cluster3", 
        "yellow",ifelse(annotations["category"]=="Cluster4", 
        "green",ifelse(annotations["category"]=="Cluster5", 
        "blue",ifelse(annotations["category"]=="Cluster6", 
        "purple",ifelse(annotations["category"]=="Cluster7", 
           "skyblue",ifelse(annotations["category"]=="Cluster8", 
        "black","grey"))))))))
    return(colorsVector)
}

testHeatmap3<-function(logCPM, annotations) {    
    sampleColors = mapDrugToColor(annotations)
    # Assign just column annotations
    heatmap3(logCPM, margins=c(10,10), 
    	ColSideColors=sampleColors,scale="none",
    	col = colorRampPalette(c( "#0072B2","#F0E442", "#D16103"))(1024),
    	Rowv=NA,
    	Colv=NA,
        xlab = "Cell ID",
    ylab = "Marker genes",
    showColDendro = F,
  showRowDendro = F) 
    #Assign column annotations and make a custom legend for them
    heatmap3(logCPM, margins=c(10,10), ColSideColors=sampleColors, 
    	scale="none",
    	col = colorRampPalette(c( "#0072B2", "#F0E442","#D16103"))(1024),
        legendfun=function()showLegend(legend=paste0("Cluster",1:7), col=c("red", "orange", "yellow","green","blue","purple","skyblue"), cex=1),
            	Rowv=NA,
    	Colv=NA,
        xlab = "Cell ID",
    ylab = "Marker genes",
    showColDendro = F,
  showRowDendro = F)
    
    #Assign column annotations as a mini-graph instead of colors,
    #and use the built-in labeling for them
    ColSideAnn<-data.frame(Cluster=annotations[["category"]])
    heatmap3(logCPM, ColSideAnn=ColSideAnn,
        #ColSideFun=function(x)showAnn(x),
        margins=c(10,10),
        ColSideWidth=0.8,
        Rowv=NA,
    	Colv=NA,
        xlab = "Cell ID",
    	ylab = "Marker genes",
    	showColDendro = F,
  		showRowDendro = F)
}


# library(BiRewire)
library(corrplot)
library(heatmap3)

marker_TLS = c("CCL2","CCL3","CCL4","CCL5","CCL8","CCL18","CCL19","CCL21","CXCL9","CXCL10",
    "CXCL11","CXCL13","CD200","FBLN7","ICOS","SGPP2","SH2D1A","TIGIT","PDCD1")
heat_expr = SCTscaled[which(rownames(SCTscaled) %in% marker_TLS),]

anno = paste0("Cluster",as.character(meta_data$SpatialPCA_Walktrap))

myorder = order(anno)
category = anno[myorder]
gAnnotationData = data.frame("cells"=colnames(heat_expr)[myorder], category)
colnames(gAnnotationData) = c("cells","category")
gLogCpmData = heat_expr[,myorder]
genenum = dim(heat_expr)[1]
datt = matrix(0,genenum,7)
for(i in 1:7){
    for(j in 1:genenum){
        datt[j,i] = mean(gLogCpmData[j,which(category %in% paste0("Cluster",i))])
    }
}

colnames(datt) = paste0(1:7)
rownames(datt) = rownames(heat_expr)
category = paste0("Cluster",1:7)
gAnnotationData_mean = data.frame("cells"=colnames(datt), category)
colnames(gAnnotationData_mean) = c("cells","category")
na.omit(match(marker_TLS,rownames(datt)))
datt = datt[rev(na.omit(match(marker_TLS,rownames(datt)))),]

pdf("TLS_marker_heatmap.pdf",width=6, height=6)  #heat_expr = scale(SCTscaled)
testHeatmap3(datt, gAnnotationData_mean)
dev.off()

```

Marker gene mean expression in each cluster plot.
```R

name = "FN1"
marker_cancer = name
heat_expr = SCTscaled[which(rownames(SCTscaled) %in% marker_cancer),]
names(heat_expr) = colnames(SCTscaled)
anno = paste0("Cluster",as.character(meta_data$SpatialPCA_Walktrap))
myorder = order(anno)
category = anno[myorder]

gAnnotationData = data.frame("cells"=names(heat_expr)[myorder], category)
colnames(gAnnotationData) = c("cells","category")
gLogCpmData = heat_expr[myorder]
mean_expr_ERBB2 = c()
for(i in 1:7){
        mean_expr_ERBB2[i] = mean(gLogCpmData[which(category %in% paste0("Cluster",i))])
}
datt = data.frame("Mean_expression"=mean_expr_ERBB2,"Cluster" = paste0("Cluster",1:7))
cbp1 <- c(  "plum1", "dodgerblue","mediumaquamarine",  "palegreen4", "chocolate1",
            "lightblue2","#F0E442",  "red","#CC79A7","mediumpurple","seagreen1")
pdf(paste0(marker_cancer,"_mean_expression_barplot.pdf"),width=5,height=5)
ggplot(data=datt, aes(x=Cluster, y=Mean_expression, fill=Cluster)) +
  geom_bar(stat="identity", position=position_dodge(),color="black",width=0.8)+
   # scale_fill_brewer(palette="Paired")+
  # scale_colour_manual(values=cbp1)+
  scale_fill_manual(values = cbp1)+
  labs(title=paste0("Mean expression of ",marker_cancer),x="", y = "")+
  theme_classic(base_size=18)+
  theme(plot.title = element_text(size = 15),
              text = element_text(size = 15),
              #axis.title = element_text(face="bold"),
              #axis.text.x=element_text(size = 20) ,
              legend.position = "right")# +
dev.off()


name = "ERBB2"
marker_cancer = name
heat_expr = SCTscaled[which(rownames(SCTscaled) %in% marker_cancer),]
names(heat_expr) = colnames(SCTscaled)
anno = paste0("Cluster",as.character(meta_data$SpatialPCA_Walktrap))
myorder = order(anno)
category = anno[myorder]

gAnnotationData = data.frame("cells"=names(heat_expr)[myorder], category)
colnames(gAnnotationData) = c("cells","category")
gLogCpmData = heat_expr[myorder]
mean_expr_ERBB2 = c()
for(i in 1:7){
        mean_expr_ERBB2[i] = mean(gLogCpmData[which(category %in% paste0("Cluster",i))])
}
datt = data.frame("Mean_expression"=mean_expr_ERBB2,"Cluster" = paste0("Cluster",1:7))
cbp1 <- c(  "plum1", "dodgerblue","mediumaquamarine",  "palegreen4", "chocolate1",
            "lightblue2","#F0E442",  "red","#CC79A7","mediumpurple","seagreen1")
pdf(paste0(marker_cancer,"_mean_expression_barplot.pdf"),width=5,height=5)
ggplot(data=datt, aes(x=Cluster, y=Mean_expression, fill=Cluster)) +
  geom_bar(stat="identity", position=position_dodge(),color="black",width=0.8)+
   # scale_fill_brewer(palette="Paired")+
  # scale_colour_manual(values=cbp1)+
  scale_fill_manual(values = cbp1)+
  labs(title=paste0("Mean expression of ",marker_cancer),x="", y = "")+
  theme_classic(base_size=18)+
  theme(plot.title = element_text(size = 15),
              text = element_text(size = 15),
              #axis.title = element_text(face="bold"),
              #axis.text.x=element_text(size = 20) ,
              legend.position = "right")# +
dev.off()


```








