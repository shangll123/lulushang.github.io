---
layout: default
title: ggplot figures
parent: Academic
nav_order: 4
---

##### violin plot 

```
pdf(paste0("Number_independent_eQTL_vs_ratio_violinplot_EA.pdf"),width=10, height=6)
dp <- ggplot(dat, aes(x=num, y=pve)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="European American",x="Number of independent eQTLs", y = "Ratio of indep eQTLs-PVE / cis-PVE")+
  scale_fill_brewer(palette="RdBu") + 
  theme_minimal(base_size = 22)
dp
dev.off()
```

##### boxplot
```
#-----
pdf(paste0("Number_independent_eQTL_vs_ratio_boxplot_AA.pdf"),width=6, height=6)
p <- ggplot(dat, aes(x = num, y = pve )) +
        geom_boxplot(alpha = 0.8,fill = "cornflowerblue") + 
        scale_fill_manual(values=c("chocolate1","seagreen3"))+
        scale_y_continuous(name = "Ratio of indep eQTLs-PVE / cis-PVE",
                           #breaks = seq(0, 9, 1),
                           limits=c(0, 1)) +
        scale_x_discrete(name = "Number of independent eQTLs") +
        ggtitle("African American") +
        theme_minimal() +
        theme(plot.title = element_text(size = 18),
              text = element_text(size = 18),
              #axis.title = element_text(face="bold"),
              axis.text.x=element_text(size = 18) ,
              legend.position = "bottom")# +
        #facet_grid(. ~ Method)
p
dev.off()

#-----

dat = data.frame(Method, ranking, Tissue, Tissue_types)
pdf(paste0(disease,"_grid_A1.pdf"),width=5, height=5)
p10 <- ggplot(dat, aes(x = Tissue_types, y = ranking, fill = Tissue_types)) +
        geom_boxplot(alpha = 0.7) +scale_fill_manual(values=c("skyblue2","seagreen3"))+
        scale_y_continuous(name = "Rank of Tissues",
                           breaks = seq(0, 40, 5),
                           limits=c(0, 40)) +
        scale_x_discrete(name = " ") +
        ggtitle(paste0(disease_name[i])) +
        theme_bw() +
        theme(plot.title = element_text(size = 20,  face = "bold"),
              text = element_text(size = 20),
              axis.title = element_text(face="bold"),
              axis.text.x=element_text(size = 15) ,
              legend.position = "bottom") +
        facet_grid(. ~ Method)
p10
dev.off()


```
##### histogram

```

pdf(paste0("Histogram_",num_variable[k],".pdf"),width=8, height=8)
ggplot(df, aes(x=df[,k])) +
  geom_histogram(bins = 150,color="darkblue", fill="lightblue")+
  geom_vline(aes(xintercept=median(df[,k])),linetype="dashed")+ 
  theme_minimal(base_size = 22)+
  geom_density(alpha=0.6)+
   labs(title=paste0(num_variable[k]),x=paste0(num_variable[k]), y = "Count")
dev.off()
```
##### density plot
```

pdf("density_eqtl_EU_200kb.pdf")
  ggplot(eqtl_table_EU2, aes(x = eqtl_table_EU_rs_dist_tss)) + 
  geom_density(aes(fill = eQTL_order), alpha = 0.4) + 
  geom_vline(data = mu, aes(xintercept = grp.median, color = eQTL_order), linetype = "dashed") + 
  labs(title="EU population, distance to tss within 200kb",x ="Distance to TSS(kb)") 
dev.off()


```
##### barplot
```
pdf("Overlap_3methods_mouseOB.pdf",width=14, height=8)
ggplot(data=dat, aes(x=topcount, y=count, fill=method)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  labs(title="Compare with spatialDE",x="Number of Top Genes selected", y = "Overlapped genes")+
  theme_minimal(base_size = 22)+
  theme(legend.position="bottom") 
dev.off()
```

##### pie plot

```

df = data.frame("State" = statenames,"Percentage" = statenum)
library(ggplot2)
pdf("Pie_state.pdf",width=8, height=8)
pie = ggplot(df, aes(x="", y=Percentage, fill=State)) + geom_bar(stat="identity", width=1)
pie = pie + coord_polar("y", start=0) + geom_text(aes(label = paste0(round(Percentage*100), "%")), position = position_stack(vjust = 0.5))
pie = pie + scale_fill_brewer(palette = "Set3") 
pie = pie + labs(x = NULL, y = NULL, fill = NULL, title = "State")
pie = pie + theme_classic() + theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0.5, color = "#666666"))
pie
dev.off()         
```

##### scatter plot
```
pdf(paste0("update_single_tissue_simu_rank_scatter_my_h",my_h,".pdf"))
ggplot(mydata, aes(x=originalranking, y=subs)) +
  geom_point(shape=19, fill="#56B4E9", color="#56B4E9", size=3)+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="tomato1")+
  theme_bw(base_size=25)+
  labs(title=paste0("sigmal strength: ",my_h),
       x="Tissue ranks in the original dataset", y = "Tissue ranks in perturbed dataset")
dev.off()

```

##### line plot
```

#--------------------------
pdf(paste0("update_single_tissue_simu_reproducibility_my_h",my_h,".pdf"))
ggplot(topss, aes(x=tis, y=reproducibility)) +
  geom_point(shape=19, fill="hotpink1", color="hotpink1", size=3)+
  geom_line() +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="tomato1")+
  theme_bw(base_size=25)+
  ylim(0,1)+
  labs(title=paste0("Signal strength: ",my_h),
       x="Tissue ranks in the original dataset", y = "Reproducibility")
dev.off()

#--------------------------
#  y: #eQTLs, x: #PCs
PC = c(0,  5, 10, 15, 20)
nG = PC_eQTL
nPC = length(PC)
pdf("eSNP_vs_PC_AA.pdf", 4,3)
plot(0, 0, xlim=c(0, 20), ylim=c(0, 600000), type="n", xlab="Number of Genotype PCs", ylab="Number of eSNPs", main="",cex.lab=1.2)
abline(v=seq(0, 20, 5), lty=2, col="lightgrey")
abline(h=seq(0, 600000, 50000), lty=2, col="lightgrey")
points(PC, nG, type="o", col=COL[2], pch=20, lwd=2)
#legend("bottomright", legend=c("eQTLs vs PCs at 5% FDR in African American"), fill=c(COL[1], "grey"), bg="white", cex=1.2)
dev.off()


```

##### heatmap
```
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

```

