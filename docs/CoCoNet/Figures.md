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



