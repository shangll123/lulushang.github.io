---
layout: default
title: Breast Tumor Analysis
nav_order: 5
has_children: false
parent: SpatialPCA
permalink: /docs/Projects/SpatialPCA/HER2ST
---

##### Load package
```R
library(SpatialPCA)
```

##### Following data processing codes from [https://github.com/almaan/her2st](https://github.com/almaan/her2st):

ST breast tumor data are downloaded from [https://github.com/almaan/her2st](https://github.com/almaan/her2st). We also saved the raw data that we used in our examples in RData format, which can be downloaded from [here](https://drive.google.com/drive/folders/1mkXV3kQKqwxk42SW4Rb263FgFj2K8HhT?usp=sharing).

```R
#   library(Seurat)
  # library(data.table)
  # library(ggplot2)
  # library(plotly)
  # library(STutility)
  # library(zeallot)
  # library(openxlsx)
  # meta_data <- read.csv("../res/ST-cluster/sample.csv",header=TRUE, stringsAsFactors=FALSE,sep=",")
# meta_data$patient_id = c()
# for(i in 1:dim(meta_data)[1]){
#   meta_data$patient_id[i] = paste0(meta_data$patient[i],meta_data$cluster[i] )
# }
# rownames(meta_data) = meta_data$patient_id
# samples <- list.files(pattern = ".tsv", path = "../data/ST-cnts/", full.names = T)
# names(samples) <- substr(do.call(rbind, strsplit(samples, split = "/"))[, 5], start = 1, stop = 2)
# imgs <- list.files(path = "../data/ST-imgs/", recursive = T, full.names = T, pattern = ".jpg")
# names(imgs) <- do.call(rbind, strsplit(imgs, split = "/"))[, 6]
# ids <- names(samples)
# infoTable <- data.frame(samples, imgs = imgs[ids], ids, patient_id = substr(x = ids, start = 1, stop = 1), stringsAsFactors = FALSE)

# tmp = meta_data[match(infoTable$ids, meta_data$patient_id),]

# infoTable <- cbind(infoTable, tmp)
# infoTable[, -c(11:28)]

# #Subset infoTable to include specified datasets.
# infoTable$spotfiles <- list.files(path = "../data/ST-spotfiles", full.names = T)[1:36]
# head(infoTable)

# ## Load data
# #Load all patient datasets and merge into one Seurat object per patient. Each gene has to bre present in at least 20 spots per sample and each spot has to have at least 300 unique features (genes).
# seu.list <- lapply(unique(infoTable$patient_id), function(s) {
#     InputFromTable(infotable = subset(infoTable, patient_id == s), 
#                       min.gene.spots = 20,
#                       min.spot.feature.count = 300,
#                       platform = "1k")
# }
# )

# # remove ring genes
# seu.list <- lapply(seu.list, function(seu) {
#   subset(seu, features = setdiff(rownames(seu@assays$RNA@counts), ring.genes))
# })


# #Calculate some QC metrics
# total.qc <- do.call(rbind, lapply(seu.list, function(se) {
#   data.frame(total_UMIs = sum(se@assays$RNA@counts), nSpots = ncol(se))
# }))


# # QC

# qcMat <- do.call(rbind, lapply(1:length(seu.list), function(i) {
#     seu <- seu.list[[i]]
#     do.call(rbind, lapply(unique(seu[["ids", drop = T]]), function(id) {
#         repMat <- seu@assays$RNA@counts[, seu[["ids", drop = T]] == id]
#         nUMI <- Matrix::colSums(repMat)
#         nGene <- apply(repMat, 2, function(x) sum(x > 0))
#         data.frame(sample = id, 
#                    avg.nUMI = round(mean(nUMI)),
#                    median.nUMI = median(nUMI),
#                    max.nUMI = max(nUMI),
#                    min.nUMI = min(nUMI),
#                    avg.nGene = round(mean(nGene)),
#                    median.nGene = median(nGene),
#                    min.nGene = min(nGene),
#                    max.nGene = max(nGene),
#                    nSpots = ncol(repMat))
#     }))
# }))

# qcMat

# ```

# Prepare count matrix and normalized matrix. We used H1 sample in the manuscript. We took the raw count matrix and location matrix as input in our paper. 
# ```R


# # default variable.features.rv.th is 1.3, original https://github.com/almaan/her2st used 1.1
# seu.list.1.3 = seu.list
# seu.list.1.3 <- lapply(seu.list.1.3, function(seu) {
#   SCTransform(seu, 
#               vars.to.regress = c("ids"), 
#               return.only.var.genes = FALSE, 
#               variable.features.n = NULL, 
#               variable.features.rv.th = 1.3)
# })

# for(num in 1:length(seu.list.1.3)){
#   print(num)
#   seu.list.single = seu.list.1.3[[num]]
#   save(seu.list.single, file = paste0("~/her2st/data/seu.list.single",num,".rv1.3.RData"))
# }

# her2stdatanum = 0
# for(num in 1:8){
#   load(paste0("~/her2st/data/seu.list.single",num,".rv1.3.RData"))
#   count = 0
#   for(id in unique(seu.list.single[["ids", drop = T]])){
#     count = count + 1
#     her2stdatanum = her2stdatanum + 1
#     print(her2stdatanum)
#     SCTcounts = seu.list.single@assays$SCT@counts[, seu.list.single[["ids", drop = T]] == id]
#     SCTscaled = seu.list.single@assays$SCT@scale.data[, seu.list.single[["ids", drop = T]] == id]
#     ind = which(seu.list.single@tools$Staffli@meta.data$sample == count)
#     metadata = seu.list.single@tools$Staffli@meta.data[ind,]
#     print(dim(metadata))
#     save(SCTcounts, SCTscaled, metadata, file = paste0("~/her2st/data/her2stdata",her2stdatanum,".rv1.3.RData"))
#   }
# }
```

###### Load data
Or you can directly download the data from above mentioned google drive.
```R
load("ST_data.RData") 
print(dim(rawcount)) # The count matrix
print(dim(location)) # The location matrix

```

##### Run SpatialPCA
```R
library(SPARK)
library(Seurat)
library(peakRAM)
library(ggplot2)
library(mclust) # ARI
library(aricode)# NMI
# location matrix: n x 2, count matrix: g x n.
# here n is spot number, g is gene number.
# here the column names of sp_count and rownames of location should be matched
ST = CreateSpatialPCAObject(counts=rawcount, location=location, project = "SpatialPCA",gene.type="spatial",sparkversion="spark", gene.number=3000,customGenelist=NULL,min.loctions = 20, min.features=20)
ST = SpatialPCA_buildKernel(ST, kerneltype="gaussian", bandwidthtype="SJ")
ST = SpatialPCA_EstimateLoading(ST,fast=FALSE,SpatialPCnum=20)
ST = SpatialPCA_SpatialPCs(ST, fast=FALSE)

## The bandwidth is:  0.048384007565889  
#   Elapsed_Time_sec Total_RAM_Used_MiB Peak_RAM_Used_MiB
# 1           55.104                -11             322.7
# > T
# Time difference of 55.10455 secs

# collect result
ind_na = which(metadataST$truth2020=="undetermined")
SpatialPCA_result = list()
SpatialPCA_result$SpatialPCs  = ST@SpatialPCs
SpatialPCA_result$normalized_expr  = ST@normalized_expr
SpatialPCA_result$location = ST@location
pred_cluster= walktrap_clustering(7, ST@SpatialPCs,round(sqrt(dim(ST@location)[1])))
SpatialPCA_result$clusterlabel = pred_cluster
SpatialPCA_result$clusterlabel_refine=refine_cluster_10x(pred_cluster,SpatialPCA_result$location,shape="square")
SpatialPCA_result$truth = metadataST$truth2020
SpatialPCA_result$ARI_original = adjustedRandIndex(SpatialPCA_result$clusterlabel[-ind_na],SpatialPCA_result$truth[-ind_na])
SpatialPCA_result$ARI = adjustedRandIndex(SpatialPCA_result$clusterlabel_refine[-ind_na],SpatialPCA_result$truth[-ind_na])
SpatialPCA_result$NMI = NMI(as.factor(SpatialPCA_result$clusterlabel_refine[-ind_na]),as.factor(SpatialPCA_result$truth[-ind_na]))
SpatialPCA_result$PAS = fx_PAS(SpatialPCA_result$clusterlabel,SpatialPCA_result$location)
SpatialPCA_result$CHAOS = fx_CHAOS(SpatialPCA_result$clusterlabel,SpatialPCA_result$location)
metadataST = data.frame("SpatialPCA" = SpatialPCA_result$clusterlabel_refine)
LISI <- compute_lisi(SpatialPCA_result$location, metadataST, c('SpatialPCA'))
SpatialPCA_result$LISI = LISI 

save(SpatialPCA_result, file = "ST_SpatialPCA_result.RData")

SpatialPCA_result$ARI
# > SpatialPCA_result$ARI
# [1] 0.4312356


# visualize the result
cbp_spatialpca <- c(  "mediumaquamarine", "lightblue2","#F0E442",  "plum1","chocolate1","dodgerblue","palegreen4","red","#CC79A7","mediumpurple","seagreen1")
pdf("ST_SpatialPCA.pdf")
clusterlabel = SpatialPCA_result$clusterlabel_refine
loc1 = ST@location[,1]
loc2 = ST@location[,2]
datt = data.frame(clusterlabel, loc1, loc2)
p = ggplot(datt, aes(x = loc1, y = loc2, color = clusterlabel)) +
            geom_point( alpha = 1,size=7) +
            scale_color_manual(values = cbp_spatialpca)+
            theme_void()+
            theme(plot.title = element_text(size = 10,  face = "bold"),
              text = element_text(size = 10),
              legend.position = "bottom")
print(p)
dev.off()

```

##### Compare Pseudo R2 of latent components in dimension reduction methods
ARI comparison figures can be obtained in a similar fasion.
```R
library(nnet)
library( DescTools  )

# Obtain PCA and NMF low dimension components
PCA_spatialgene = get_PCA(SpatialPCA_result$normalized_expr,20)
NMF_spatialgene = get_NMF(rawcount[match(rownames(SpatialPCA_result$normalized_expr),rownames(rawcount)),],20)
# read in latent embeddings from SpaGCN
SpaGCN_embed = read.csv("~/SpaGCN_her2st_embed.csv",header=F)


SpatialPCA_R2=c()
PCA_R2=c()
NMF_R2=c()
SpaGCN_R2=c()
# need to remove undetermined spots
for(i in 1:20){
	if(i==1){
		fit <- multinom(SpatialPCA_result$truth[-ind_na] ~ SpatialPCA_result$SpatialPCs[1:i,-ind_na], maxit = 1000, MaxNWts = 2000,model = TRUE)
		SpatialPCA_R2[i] = PseudoR2(fit,c("McFaddenAdj"))
		fit <- multinom(SpatialPCA_result$truth[-ind_na] ~ PCA_spatialgene[1:i,-ind_na], maxit = 1000, MaxNWts = 2000,model = TRUE)
		PCA_R2[i] = PseudoR2(fit,c("McFaddenAdj"))
		fit <- multinom(SpatialPCA_result$truth[-ind_na]~ SpaGCN_embed[-ind_na,1:i], maxit = 1000, MaxNWts = 2000,model = TRUE)
		SpaGCN_R2[i] = PseudoR2(fit,c("McFaddenAdj"))
		fit <- multinom(SpatialPCA_result$truth[-ind_na]~ NMF_spatialgene[1:i,-ind_na], maxit = 1000, MaxNWts = 2000,model = TRUE)
		NMF_R2[i] = PseudoR2(fit,c("McFaddenAdj"))

	}else{
		fit <- multinom(SpatialPCA_result$truth[-ind_na] ~ as.matrix(as.data.frame(t(SpatialPCA_result$SpatialPCs[1:i,-ind_na]))), maxit = 1000, MaxNWts = 2000,model = TRUE)
		SpatialPCA_R2[i] = PseudoR2(fit,c("McFaddenAdj"))
		fit <- multinom(SpatialPCA_result$truth[-ind_na] ~ as.matrix(as.data.frame(t(PCA_spatialgene[1:i,-ind_na]))), maxit = 1000, MaxNWts = 2000,model = TRUE)
		PCA_R2[i] = PseudoR2(fit,c("McFaddenAdj"))
		fit <- multinom(SpatialPCA_result$truth[-ind_na] ~ as.matrix(as.data.frame(SpaGCN_embed[-ind_na,1:i])), maxit = 1000, MaxNWts = 2000,model = TRUE)
		SpaGCN_R2[i] = PseudoR2(fit,c("McFaddenAdj"))
		fit <- multinom(SpatialPCA_result$truth[-ind_na] ~ as.matrix(as.data.frame(t(NMF_spatialgene[1:i,-ind_na]))), maxit = 1000, MaxNWts = 2000,model = TRUE)
		NMF_R2[i] = PseudoR2(fit,c("McFaddenAdj"))
	}
}

fit <- multinom(SpatialPCA_result$truth[-ind_na] ~ metadataST$BayesSpace_hvggene[-ind_na], maxit = 1000, MaxNWts = 2000,model = TRUE)
BayesSpace_R2 = PseudoR2(fit,c("McFaddenAdj"))

BayesSpace_R2 = rep(BayesSpace_R2, 20)
PseudoR2=c(SpatialPCA_R2, PCA_R2, NMF_R2,SpaGCN_R2,BayesSpace_R2)
method=c(rep("SpatialPCA",20),rep("PCA",20),rep("NMF",20),rep("SpaGCN",20),rep("BayesSpace",20))
PCnum=c(1:20,1:20,1:20,1:20,1:20)
method=factor(method, levels=c("SpatialPCA","SpaGCN","BayesSpace","PCA","NMF"))
dat=data.frame(PseudoR2, method, PCnum)
pdf(paste0("ST_addingPCs_PseudoR2_jan16.pdf"),width=10,height=5)
p<- ggplot(dat, aes(x=PCnum, y=PseudoR2, color=method,group=method)) + 
   theme_bw(base_size = 22)+
   geom_line(size=1)+
  # ylim(0,1)+
  geom_point( size=2, color="black")+
  labs(title=paste0("Pseudo R-square (Tumor data)"),x="Low dimensional components number", y = "Pseudo R2")+
   theme(plot.title = element_text(size = 20),
              text = element_text(size = 20),
              legend.position = "right")# +
print(p)
dev.off()

```

##### Trajectory analysis
In ST data, we focused on trajectory inference on tumor and tumor adjacent regions to investigate how these locations are connected to one another and underlie tumorigenesis. Based on the inferred pseudo-time values, we connected neighboring locations on the tissue to construct trajectories.
```R

library(slingshot)
# focus on tumor and its surrounding region
tumor_ind = which(clusterlabel_refine %in% c(5,4,6))
sim<- SingleCellExperiment(assays = rawcount[,tumor_ind])
reducedDims(sim) <- SimpleList(DRM = t(ST@SpatialPCs[,tumor_ind]))
colData(sim)$Walktrap <- factor(clusterlabel_refine[tumor_ind])    
sim  <-slingshot(sim, clusterLabels = 'Walktrap', reducedDim = 'DRM',start.clus="5" ) 
# in this data we set tumor region as start cluster
summary(sim@colData@listData)

pseudotime_traj1 = sim@colData@listData$slingPseudotime_1
gridnum = 10
cbp <- c(  "plum1", "chocolate1","dodgerblue","black")
p_traj1 = plot_trajectory(pseudotime_traj1, ST@location[tumor_ind,],clusterlabel_refine[tumor_ind],gridnum,cbp,pointsize=5 ,arrowlength=0.2,arrowsize=1,textsize=15 )
p_traj1$Arrowoverlay1
p_traj1$Pseudotime

```

##### High resolution map reconstruction
```R
STsimu_high_ST = SpatialPCA_highresolution(ST)
cluster_SpatialPCA_high = walktrap_clustering(7, latent_dat=STsimu_high_ST@highPCs,76)
color_in=c(  "plum1", "palegreen4","mediumaquamarine",  "chocolate1","#F0E442","dodgerblue","lightblue2")
title_in="SpatialPCA High resolution"
plot_cluster(STsimu_high_ST@highPos, as.character(cluster_SpatialPCA_high), pointsize=2,text_size=20 ,title_in,color_in,legend="bottom")

# By default the new unmeasured locations are generated based on the existing spots. We also generated custom locations and predicted the high resolution map.
loc_round = round(ST@location/100) # because we are using pixel coordinates provided by the original data paper, the scale of location coordinates is large, we divided by 100 for easily generating spots that cover the whole tissue slice.

x_in = c()
y_in = c()
for(y_coor in unique(loc_round[,2])){
    ind = which(loc_round[,2]==y_coor)
    x_ind = loc_round[ind,1]
    x_new = seq(from=min(x_ind),to=max(x_ind),by=1)
    y_new = rep(y_coor, length(x_new))
    x_in = c(x_in,x_new)
    y_in = c(y_in, y_new)
}
loc_round=data.frame(x_in,y_in)
x_in = c()
y_in = c()
for(x_coor in unique(loc_round[,1])){
    ind = which(loc_round[,1]==x_coor)
    y_ind = loc_round[ind,2]
    y_new = seq(from=min(y_ind),to=max(y_ind),by=1)
    x_new = rep(x_coor, length(y_new))
    y_in = c(y_in,y_new)
    x_in = c(x_in, x_new)
}
dat_sample = data.frame(x_in,y_in)*100

STsimu_high_ST_random = SpatialPCA_highresolution(ST,newloc = dat_sample)
cluster_SpatialPCA_high_custom = walktrap_clustering(8, latent_dat=STsimu_high_ST_random@highPCs,round(sqrt(dim(STsimu_high_ST_random@highPCs)[2])))
color_in=c(  "dodgerblue", "mediumaquamarine","#CC79A7",  "chocolate1","lightblue2","palegreen4","plum1","#F0E442","#CC79A7","mediumpurple","seagreen1")
title_in="SpatialPCA High resolution"
plot_cluster(STsimu_high_ST_random@highPos, as.character(cluster_SpatialPCA_high_custom), pointsize=1.5,text_size=20 ,title_in,color_in,legend="bottom")

```

##### Differential gene analysis with Seurat
```R
Seu <- CreateSeuratObject(counts = rawcount, project = "ST", min.cells = 20, min.features = 20)
Seu = SCTransform(Seu, return.only.var.genes = FALSE, variable.features.n = NULL,  variable.features.rv.th = 1.3)
Seu <- RunPCA(Seu, features = VariableFeatures(object = Seu))
Seu <- FindNeighbors(Seu, dims = 1:10)
Seu <- FindClusters(Seu, resolution = 0.5)

Idents(Seu) = paste0("cluster",SpatialPCA_result$clusterlabel_refine)
DE_gene = list()
DEgene_spatialPCA=c()
each_num = c()
for(cluster in 1:7){
	print(cluster)
	DE_gene[[cluster]] = FindMarkers(Seu, ident.1 = paste0("cluster",cluster), ident.2 =NULL, test.use = "MAST")
	each_num[cluster] = dim(DE_gene[[cluster]])[1]
	DEgene_spatialPCA = c(DEgene_spatialPCA, rownames(DE_gene[[cluster]]))
}
DEgene_spatialPCA=unique(DEgene_spatialPCA)
length(DEgene_spatialPCA)

# > length(DEgene_spatialPCA)
# [1] 940

each_num
> each_num
[1] 143 104 321 131 389  81 369
```

##### Marker gene mean expression box plot
```R
SCTscaled = Seu@assays$SCT@scale.data
name = "DDX54"
name = "TRAF2"

marker_cancer = name
heat_expr = SCTscaled[which(rownames(SCTscaled) %in% marker_cancer),]
names(heat_expr) = colnames(SCTscaled)
anno = paste0("Cluster",as.character(metadataST$SpatialPCA))
myorder = order(anno)
category = anno[myorder]

gAnnotationData = data.frame("cells"=names(heat_expr)[myorder], category)
colnames(gAnnotationData) = c("cells","category")
gLogCpmData = heat_expr[myorder]
mean_expr_ERBB2 = c()
for(i in 1:7){
        mean_expr_ERBB2[i] = mean(gLogCpmData[which(category %in% unique(paste0("Cluster",as.character(metadataST$SpatialPCA)))[i])])
}

truths=unique(paste0(as.character(metadataST$SpatialPCA)))
# truths = factor(truths, levels=c("connective tissue","cancer in situ","invasive cancer","immune infiltrate","breast glands","adipose tissue","undetermined"),order=T)
datt = data.frame("Mean_expression"=mean_expr_ERBB2,"Cluster" = truths)
cbp1=c(  "mediumaquamarine", "lightblue2","#F0E442",  "plum1","chocolate1","dodgerblue","palegreen4","red","#CC79A7","mediumpurple","seagreen1")
pdf(paste0(marker_cancer,"_mean_expression_barplot_SpatialPCA.pdf"),width=5,height=5)
ggplot(data=datt, aes(x=Cluster, y=Mean_expression, fill=Cluster)) +
  geom_bar(stat="identity", position=position_dodge(),color="black",width=0.8)+
   # scale_fill_brewer(palette="Paired")+
  # scale_colour_manual(values=cbp1)+
  scale_fill_manual(values = cbp1)+
  labs(title=paste0("SpatialPCA: Mean expression of ",marker_cancer),x="", y = "")+
  theme_bw(base_size=18)+
  ylim(-0.25,0.6)+
  theme(plot.title = element_text(size = 15),
              text = element_text(size = 15),
              #axis.title = element_text(face="bold"),
              #axis.text.x=element_text(size = 20) ,
              legend.position = "right")# +
dev.off()


```

##### TLS region identification
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

# use all genes normalized expression to make marker gene heatmap

seu <- CreateSeuratObject(counts = rawcount, project = "forheatmap", min.cells = 0, min.features = 0)
seu=   SCTransform(seu,          
              return.only.var.genes = FALSE, 
              variable.features.n = NULL, 
              variable.features.rv.th = 1.3)     
SCTscaled = seu@assays$SCT@scale.data

marker_TLS = c("CCL2","CCL3","CCL4","CCL5","CCL8","CCL18","CCL19","CCL21","CXCL9","CXCL10",
    "CXCL11","CXCL13","CD200","FBLN7","ICOS","SGPP2","SH2D1A","TIGIT","PDCD1")
heat_expr = SCTscaled[which(rownames(SCTscaled) %in% marker_TLS),]

anno = paste0("Cluster",as.character(metadataST$SpatialPCA))

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

##### Deconvolution of cell types
We performed two cell type deconvolution analysis, one with immune cells, and another one with immune cells + cancer cells.
```R

#=====================================================================================
# First Deconvolution of cell types, with immune cells
#=====================================================================================

# devtools::install_github("dmcable/RCTD", build_vignettes = TRUE)
library(RCTD)
library(Matrix)
# GSE114725
# - Azizi, Elham, Ambrose J. Carr, George Plitas, Andrew E. Cornish, Catherine Konopacki, Sandhya Prabhakaran, Juozas Nainys, et al. “[Single-Cell Map of Diverse Immune Phenotypes in the Breast Tumor Microenvironment](https://doi.org/10.1016/j.cell.2018.05.060).” Cell, June 2018. - scRNA-seq of immune cells in BRCA - continuous activation of T cells, no macrophage polarization. inDrop and 10X platforms. 47,016 CD45+ cells from 8 primary breast carcinomas. 83 clusters, tested by cross-validation. [Data 1](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114727), [Data 2](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114725), [Data 3](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114724)
GSE114725 = read.csv("~/raw_corrected.csv",header=T)
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
celltype = as.character(results$results_df$first_type)
metadataST$celltype = celltype

#=====================================================================================
# Second Deconvolution of cell types, with immune cells + cancer cells
#=====================================================================================

# Wu et al. EMBO 2020
# Stromal cell diversity associated with immune evasion in human triple-negative breast cancer
# https://www.embopress.org/doi/epdf/10.15252/embj.2019104063
# found at https://singlecell.broadinstitute.org/single_cell/study/SCP1106/stromal-cell-diversity-associated-with-immune-evasion-in-human-triple-negative-breast-cancer#study-download

library(Matrix)
library(readr)
# Read in `matrix.mtx`
counts <- readMM("./counts_matrix/matrix.mtx.gz")
# Read in `genes.tsv`
genes <- read_tsv("./counts_matrix/features.tsv.gz", col_names = FALSE)
gene_ids <- genes$X1
# Read in `barcodes.tsv`
cell_ids <- read_tsv("./counts_matrix/barcodes.tsv.gz", col_names = FALSE)$X1
rownames(counts) <- gene_ids
colnames(counts) <- cell_ids
meta_ref = read.csv("Wu_EMBO_metadata.csv")
meta_ref_singlecell = meta_ref[-1,]

 library(RCTD)
# reference single cell data
meta_data <- meta_ref_singlecell # load in meta_data (barcodes, clusters, and nUMI)
cell_types <- meta_data$celltype_final; names(cell_types) <- as.character(meta_data$NAME) # create cell_types named list
cell_types <- as.factor(cell_types) # convert to factor data type
nUMI <- as.integer(meta_data$nCount_RNA); names(nUMI) <- as.character(meta_data$NAME)  # create nUMI named list
### Create the Reference object
reference <- Reference(counts, cell_types, nUMI)
counts <- rawcount
coords <- as.data.frame(location)
nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI
### Create SpatialRNA object
puck <- SpatialRNA(coords, counts, nUMI)
# run RCTD
myRCTD <- create.RCTD(puck, reference, max_cores = 10)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
results <- myRCTD@results
celltype = as.character(results$results_df$first_type)
metadataST$celltype = celltype


#=====================================================================================
# visualize cell type proportion distribution in spatial domains
#=====================================================================================
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(ggthemes)

cbp2 = c( "#C4961A", "#F4EDCA", "deepskyblue1", "cadetblue3","lightblue1","plum1","chartreuse3","tomato","#C3D7A4",  "#52854C","#D16103","#4E84C4","#FFDB6D")

preparedata = function(percentage){
rownames(percentage) = paste0("Cluster",1:7)
celltype =  names(table( metadataST$celltype))
colnames(percentage) = celltype
rownames(percentage) = paste0("Cluster",1:7)
percentage_vec = c(percentage)
cluster_vec = c(rep(c(paste0("Cluster",1:7)),length(celltype)))
CellType = c(rep(celltype,each=7))
datt = data.frame(cluster_vec, percentage_vec,CellType)
datt$cluster_vec = factor(cluster_vec, level=paste0("Cluster",1:7))
return(datt)
}

makefigure = function(datt){
p=ggplot(datt, aes(y = percentage_vec,
             x = factor(cluster_vec ), fill = CellType)) +        ## global aes
  scale_fill_manual(values=cbp2)+
  geom_bar(position="stack", stat="identity",width=0.7,color="black") +
  theme_bw()+xlab("")+ylab("")+
  theme(plot.title = element_text(size = 20),
              text = element_text(size = 20),
              #axis.title = element_text(face="bold"),
              #axis.text.x=element_text(size = 12,angle = 60,hjust = 1) ,
              #axis.text.x=element_blank(),
              legend.position = "right")# +
return(p)
}

metadataST$celltype = as.factor(metadataST$celltype)

method="SpatialPCA"
percentage = matrix(0,7,13)
for(k in 1:7){
metadata_sub = metadataST[which(metadataST$SpatialPCA_walktrap==k ),]
match_type = metadata_sub$celltype
percentage[k,] = round(unlist(table(match_type))/dim(metadata_sub)[1]*100,2)
}
datt=preparedata(percentage)
pdf(paste0("Stack_barplot_",method,".pdf"),width=8,height=5)
makefigure(datt)+ggtitle(paste0(method))
percentage = matrix(0,7,13)
for(k in 1:7){
metadata_sub = metadataST[which(metadataST$SpatialPCA_walktrap==k ),]
match_type = metadata_sub$celltype
percentage[k,] = round(unlist(table(match_type))/dim(metadataST)[1]*100,2)
}
datt=preparedata(percentage)
makefigure(datt)+ggtitle(paste0(method))
dev.off()

```

##### Pseudotime associated gene GSEA
```R
library(parallel)
library(MASS)

count_use=rawcount
pseudotime = sim_SpatialPCA@colData@listData[[3]]

pseudotime_gene = function(ind){
expr = count_use[ind,]
dat = data.frame(expr,pseudotime)
colnames(dat) = c("expr","pseudotime")
dat = na.omit(dat)
res = try(summary(m1 <- glm.nb(expr ~ pseudotime, data = dat)))
if(isTRUE(class(res)=="try-error")) { next } else { 
		p_value_traj =  res$coefficients[2,4] 
		gene = rownames(count_use)[ind]
} 
if(p_value_traj<0.05/length(expr)){
	return(gene)
}else {
	return(NA)
}
}

results = unlist(mclapply(1:dim(count_use)[1], pseudotime_gene, mc.cores = 50))
errorid = which(results=="Error in FUN(X[[i]], ...) : no loop for break/next, jumping to top level\n")
results=results[-errorid]
ST_genelist = na.omit(results)

# > length(ST_genelist)
# [1] 1713

library(gprofiler2)
library(forcats)

upload_GMT_file(gmtfile = "~/h.c24557all.v7.4.symbols.filter.update.gmt")
# gp__eeZM_IE8w_Kp8

strsplit_func = function(input){
  strsplit(input, split = "_")[[1]][1]
}



datt = list()
for(num in 1:1){
  print(num)
topGOnum = 10
bg = rownames(rawcount)
gostres <- gost(query = ST_genelist, 
                organism = "gp__eeZM_IE8w_Kp8", ordered_query = FALSE, 
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
for(num in 1:1){
  print(num)
p[[num]]=ggplot(data=datt[[num]], aes(x=fct_reorder(tolower(term_id),log10p), y=log10p, fill=Source,label = ifelse(significant ==TRUE, "*",""))) +
  geom_bar(stat="identity", position=position_dodge(),color="black",width=0.8)+
  #scale_fill_continuous(low='#F0E442', high='red', guide=guide_colorbar(reverse=TRUE)) + 
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7"))+
  labs(title=paste0("Pseudotime trajectory",num),x="Biological terms", y = "-log10(p value)")+
  coord_flip()+
  theme_classic()+
  #geom_text(vjust = 0.5) +
  geom_text(vjust = 1, nudge_y = 0.5)+
  #ylim(0,1)+
    theme(plot.title = element_text(size = 30,color="black",face="bold"),
              text = element_text(size = 30,color="black",face="bold"),
              #axis.title = element_text(size = 25,color="black",face="bold"),
              axis.text.x=element_text(size = 30,color="black",face="bold") ,
              legend.position = "right")# +
}

pdf(paste0("GSEA_SpatialPCA_msigdb_HER2ST_pseudotime_gene.pdf"),width=20,height=5)
for(num in 1:1){
  print(p[[num]])
}
dev.off()



```


##### Region specific gene GSEA
```R

upload_GMT_file(gmtfile = "~/h.c24557all.v7.4.symbols.filter.update.gmt")
# gp__EG79_4imK_ul4


# make figures
datt = list()
for(num in 1:7){
  print(num)
topGOnum = 10
bg = rownames(rawcount)
gostres <- gost(query = rownames(DE_gene[[num]]), 
                organism = "gp__EG79_4imK_ul4", ordered_query = FALSE, 
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

spatial_domain = c("Fibrous tissue near normal glands","Fat tissue","Immune region","Tumor surrounding region",
"Tumor region","Fibrous tissue near tumor",
  "Normal glands")


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
    theme(plot.title = element_text(size = 30,color="black",face="bold"),
              text = element_text(size = 30,color="black",face="bold"),
              #axis.title = element_text(size = 25,color="black",face="bold"),
              axis.text.x=element_text(size = 30,color="black",face="bold") ,
              legend.position = "right")# +
}

pdf(paste0("GSEA_region_specific_ST.pdf"),width=20,height=5)
for(num in 1:7){
  print(p[[num]])
}
dev.off()
```









