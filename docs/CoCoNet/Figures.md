---
layout: default
title: Figure
parent: CoCoNet
nav_order: 5
---

##### I used the theme_Publication, scale_fill_Publication, and scale_colour_Publication, the three functions in [ggplot theme for publication ready Plots](https://rpubs.com/Koundy/71792).

```R
library(ggplot2)
library(gridExtra)
theme_Publication <- function(base_size=14, base_family="helvetica") {  
# please see the above link for the full functions
...
}
scale_fill_Publication <- function(...){
...
}
scale_colour_Publication <- function(...){
...
}
```

##### Plot ranking of likelihoods
```R
#----------------------
# Tissue likelihood ranking 
#----------------------

library(cowplot)
library(ggplot2)
library(gridExtra) 
library(ggpubr)      

trait_current = c(1:8)
myarr=c( "SCZ","BIP", "BIPSCZ" ,"Alzheimer" , "PBC", "CD" , "UC","IBD"   )	
disease_name = c("Schizophrenia","Bipolar Disorder","Bipolar Disorder/Schizophrenia","Alzheimer's Disease","Primary biliary cholangitis","Crohn's Disease","Ulcerative colitis","Inflammatory bowel disease")


dat = list()
p = list()
count = 0
for(j in trait_current){
count = count + 1
lik_expon = c()
	for(tissue in 1:38){
	load(paste0("~path/result/Tissue_trait_",j,"_tissue_",tissue,".RData"))
		re_expon=res1
		lik_expon[tissue] = re_expon$value
	}

trait = colnames(sigma_sona)[j]
Tissues = as.factor(tissue_name)
Tissues_name = tissue_name
index = c(1:38)
tissue_color = rep(NA,38)
tissue_color[7:9] = "Brain"
tissue_color[c(13,14,31)] = "Colon & Intestine"
tissue_color[-c(7:9,13,14,31 )]="Other tissues"
lik_expon_minus_mean = lik_expon - mean(lik_expon)
ranking = rank(lik_expon_minus_mean)
mydata = data.frame( ranking, h_expon,lik_expon,Tissues_name,index,tissue_color,lik_expon_minus_mean)
mydata$lik_expon = -mydata$lik_expon # this depends on the outcome of the coconet function is likelihood or -likelihood
mydata$lik_min = mydata$lik_expon - min(mydata$lik_expon)
mydata$Trait = rep(colnames(sigma_sona)[j],dim(mydata)[1])

dat[[count]] =  mydata
p[[count]] = ggbarplot(mydata, x = "Tissues_name", y = "lik_min",
          fill = "tissue_color",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "asc",          # Sort the value in descending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Log Likelihood - minimum",
          xlab = "Tissues",
          legend.title = "GTEx Tissues",
          rotate = TRUE,
          ggtheme = theme_minimal(),
          title = disease_name[count]
          )     
}       

k=0
# Repeat codes below: 

k=k+1
tiff(paste0("Tissue_",myarr[k],".tiff"), units="in",width=10, height=10,  res=100)
grid.arrange((p[[k]] + scale_fill_Publication() +theme_Publication()),nrow=1)
dev.off()


#----------------------
# Cell likelihood ranking
#----------------------


cell_type =  c("Astrocytes","Endothelial cells","GABAergic interneurons","Microglia","Neuronal stem cells",
"Oligodendrocytes","Oligodendrocyte precursor cells","Pyramidal neurons","Granule neurons","Glutamatergic neurons")
p <- list()
count=0
dat=list()
for(j in c(1:4)){
count=count+1
lik_expon = c()
	for(tissue in 1:10){
		load(paste0("~path/result/cell_trait_",j,"_tissue_",tissue,".RData"))
		re_expon = res1
		lik_expon[tissue] = re_expon$value
	}
h_expon = h_expon
lik_expon = lik_expon
trait = colnames(sc_sigma_43)[j]
Tissues = as.factor(cell_type)
Tissues_name = cell_type
index = c(1:10)
tissue_color = rep(NA,10)
tissue_color[c(1,4,6,7)] = "Glias"
tissue_color[c(2,5)] = "Others"
tissue_color[c(3,8,9,10)] = "Neurons"
lik_expon_minus_mean = lik_expon - mean(lik_expon)
mydata = data.frame(h_expon,lik_expon,Tissues_name,index,tissue_color,lik_expon_minus_mean)
mydata$ranking = rank(mydata$lik_expon)
mydata$lik_expon = -mydata$lik_expon
mydata$lik_min = mydata$lik_expon - min(mydata$lik_expon)
mydata$Trait = rep(colnames(sigma_sona)[j],dim(mydata)[1])

dat[[count]] =  mydata
p[[count]] = ggbarplot(mydata, x = "Tissues_name", y = "lik_min",
          fill = "tissue_color",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "asc",          # Sort the value in descending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = "Log Likelihood - minimum",
          xlab = "Cell types",
          legend.title = "Cell types",
          rotate = TRUE,
          ggtheme = theme_minimal(base_size=25),
          title = disease_name[count]
          )
          
}

k=0
# Repeat codes below:

k=k+1
tiff(paste0("Cell_",myarr[k],".tiff"), units="in",width=8, height=8,  res=100)
grid.arrange((p[[k]] + scale_fill_Publication() +theme_Publication()),nrow=1)
dev.off()


```




##### Heatmap of Jaccard index for tissues / cell types

```R
library(BiRewire)
library(corrplot)
library(heatmap3)

#---------------------------------
# Tissues 
#---------------------------------

load("tissue_net.RData")
load("tissue_name.RData")
Tissue_network = tissue_net

# calculate Jaccard Index between each pair of tissues

a = matrix(0,38,38)
for(i in 1:38){
	for(j in (i+1):38){
		a[i,j] = birewire.similarity( Tissue_network[[i]],Tissue_network[[j]])
		a[j,i] = a[i,j]
	}
	print(i)
	print(a[i,])
}
colnames(a) = tissue_name
rownames(a) = tissue_name

# use heatmap3
mapDrugToColor<-function(annotations){
    colorsVector = ifelse(annotations["category"]=="Others", 
        "blue", ifelse(annotations["category"]=="Brain related", 
        "green", "red"))
    return(colorsVector)
}
testHeatmap3<-function(logCPM, annotations) {    
    sampleColors = mapDrugToColor(annotations)
    
    # Assign just column annotations
    heatmap3(logCPM, margins=c(16,16), ColSideColors=sampleColors,scale="none") 
    # Assign column annotations and make a custom legend for them
    heatmap3(logCPM, margins=c(16,16), ColSideColors=sampleColors, scale="none",
        legendfun=function()showLegend(legend=c("Others", 
        "Brain related", "Colon related"), col=c("blue", "green", "red"), cex=1))
    
    # Assign column annotations as a mini-graph instead of colors,
    # and use the built-in labeling for them
    ColSideAnn<-data.frame(Drug=annotations[["category"]])
    heatmap3(logCPM,ColSideAnn=ColSideAnn,
        ColSideFun=function(x)showAnn(x),
        ColSideWidth=0.8)
}
category = c(rep("Others",6),rep("Brain related",3),rep("Others",3),"Colon related", 
"Colon related", rep("Others",16),"Colon related","Others","Colon related",rep("Others",5))
gAnnotationData = data.frame(tissue_name, category)
gLogCpmData = a
pdf("Tissue_GTEx_heatmap_unscaled.pdf",width=7, height=7)
diag(gLogCpmData)=1    
testHeatmap3(gLogCpmData, gAnnotationData)
dev.off()

#---------------------------------
# Cells - same cell type combined
#---------------------------------
load("cell_net.RData")
Tissue_network_cell = cell_net
cell_types =  c("ASC","END","GABA","MG","NSC","ODC","OPC","exCA","exDG","exPFC")
cellnames = c("Astrocytes","Endothelial","GABAergic neurons","Microglia","Neuronal stem cells","Oligodendrocytes","Oligodendrocyte precursor cells","Pyramidal neurons","Granule neurons","Glutamatergic neurons")
category = c("Glias","Endothelial","Neurons","Glias","Neuronal Stem","Glias","Glias","Neurons","Neurons","Neurons")
a = matrix(0,10,10)
for(i in 1:10){
	for(j in (i+1):10){
		a[i,j] = birewire.similarity( Tissue_network_cell[[i]],Tissue_network_cell[[j]])
		a[j,i] = a[i,j]
	}
	print(i)
	print(a[i,])
}
colnames(a) = cellnames
rownames(a) = cellnames
diag(a)=1
mapDrugToColor<-function(annotations){
    colorsVector = ifelse(annotations["category"]=="Glias", 
        "blue", ifelse(annotations["category"]=="Neurons", 
        "green", "red"))
    return(colorsVector)
}
testHeatmap3<-function(logCPM, annotations) {    
    sampleColors = mapDrugToColor(annotations)
    # Assign just column annotations
    heatmap3(logCPM, margins=c(16,16), ColSideColors=sampleColors,scale="none") 
    # Assign column annotations and make a custom legend for them
    heatmap3(logCPM, margins=c(16,16), ColSideColors=sampleColors, ,scale="none",
        legendfun=function()showLegend(legend=c("Glias", 
        "Neurons", "Others"), col=c("blue", "green", "red"), cex=1.5))
    
    # Assign column annotations as a mini-graph instead of colors,
    # and use the built-in labeling for them
    ColSideAnn<-data.frame(Drug=annotations[["category"]])
    heatmap3(logCPM,ColSideAnn=ColSideAnn,
        ColSideFun=function(x)showAnn(x),
        ColSideWidth=0.8)
}
gAnnotationData = data.frame(cell_types, category)
gLogCpmData = a
pdf("scGTEx_heatmap_combine_unscaled.pdf",width=7, height=7)
diag(gLogCpmData)=1    
testHeatmap3(gLogCpmData, gAnnotationData)
dev.off()
```



