---
layout: default
title: Figure codes
parent: Resources
nav_order: 2
permalink: /docs/Resources/Figure
---

## violin plot 

```R
pdf(paste0("xxx.pdf"),width=10, height=6)
dp <- ggplot(dat, aes(x=num, y=pve)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(title="xxx",x="xxx", y = "xxx")+
  scale_fill_brewer(palette="RdBu") + 
  theme_minimal(base_size = 22)
dp
dev.off()
```

## boxplot
```R
#-----
pdf(paste0("xxx.pdf"),width=6, height=6)
p <- ggplot(dat, aes(x = num, y = pve )) +
        geom_boxplot(alpha = 0.8,fill = "cornflowerblue") + 
        scale_fill_manual(values=c("chocolate1","seagreen3"))+
        scale_y_continuous(name = "xxx",
                           #breaks = seq(0, 9, 1),
                           limits=c(0, 1)) +
        scale_x_discrete(name = "xxx") +
        ggtitle("xxx") +
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
pdf(paste0("xxx.pdf"),width=5, height=5)
p10 <- ggplot(dat, aes(x = Tissue_types, y = ranking, fill = Tissue_types)) +
        geom_boxplot(alpha = 0.7) +scale_fill_manual(values=c("skyblue2","seagreen3"))+
        scale_y_continuous(name = "xxx",
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
## histogram

```R

pdf(paste0("xxx.pdf"),width=8, height=8)
ggplot(df, aes(x=df[,k])) +
  geom_histogram(bins = 150,color="darkblue", fill="lightblue")+
  geom_vline(aes(xintercept=median(df[,k])),linetype="dashed")+ 
  theme_minimal(base_size = 22)+
  geom_density(alpha=0.6)+
   labs(title=paste0(xxx),x=paste0(xxx), y = "xxx")
dev.off()
```
## density plot
```R

pdf("xxx.pdf")
  ggplot(eqtl_table_EU2, aes(x = eqtl_table_EU_rs_dist_tss)) + 
  geom_density(aes(fill = eQTL_order), alpha = 0.4) + 
  geom_vline(data = mu, aes(xintercept = grp.median, color = eQTL_order), linetype = "dashed") + 
  labs(title="xxx",x ="xxx") 
dev.off()


```
## barplot
```R
pdf("xxx.pdf",width=14, height=8)
ggplot(data=dat, aes(x=topcount, y=count, fill=method)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  labs(title="xxx",x="xxx", y = "xxx")+
  theme_minimal(base_size = 22)+
  theme(legend.position="bottom") 
dev.off()
```

## pie plot

```R

df = data.frame("State" = statenames,"Percentage" = statenum)
library(ggplot2)
pdf("xxx.pdf",width=8, height=8)
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

## scatter plot
```R
pdf(paste0("xxx.pdf"))
ggplot(mydata, aes(x=originalranking, y=subs)) +
  geom_point(shape=19, fill="#56B4E9", color="#56B4E9", size=3)+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="tomato1")+
  theme_bw(base_size=25)+
  labs(title=paste0("xxx"),
       x="xxx", y = "xxx")
dev.off()

```

## line plot
```R

#--------------------------
pdf(paste0("xxx.pdf"))
ggplot(topss, aes(x=tis, y=reproducibility)) +
  geom_point(shape=19, fill="hotpink1", color="hotpink1", size=3)+
  geom_line() +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="tomato1")+
  theme_bw(base_size=25)+
  ylim(0,1)+
  labs(title=paste0("xxx"),
       x="xxx", y = "xxx")
dev.off()

#--------------------------
#  y: #eQTLs, x: #PCs
PC = c(0,  5, 10, 15, 20)
nG = PC_eQTL
nPC = length(PC)
pdf("xxx.pdf", 4,3)
plot(0, 0, xlim=c(0, 20), ylim=c(0, 600000), type="n", xlab="Number of Genotype PCs", ylab="Number of eSNPs", main="",cex.lab=1.2)
abline(v=seq(0, 20, 5), lty=2, col="lightgrey")
abline(h=seq(0, 600000, 50000), lty=2, col="lightgrey")
points(PC, nG, type="o", col=COL[2], pch=20, lwd=2)
#legend("bottomright", legend=c("eQTLs vs PCs at 5% FDR in African American"), fill=c(COL[1], "grey"), bg="white", cex=1.2)
dev.off()


```

## heatmap
```R
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

## heatmap change margin size and color
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
    heatmap3(logCPM, margins=c(10,10), ColSideColors=sampleColors,scale="none",col = colorRampPalette(c("firebrick", "yellow", "white"))(1024)) 
    # Assign column annotations and make a custom legend for them
    heatmap3(logCPM, margins=c(10,10), ColSideColors=sampleColors, scale="none",col = colorRampPalette(c("firebrick", "yellow", "white"))(1024),
        legendfun=function()showLegend(legend=c("Others", 
        "Brain related", "Colon related"), col=c("blue", "green", "red"), cex=1))
    
    # Assign column annotations as a mini-graph instead of colors,
    # and use the built-in labeling for them
    ColSideAnn<-data.frame(Drug=annotations[["category"]])
    heatmap3(logCPM,ColSideAnn=ColSideAnn,
        ColSideFun=function(x)showAnn(x),
        margins=c(10,10),
        ColSideWidth=0.8)
}
category = c(rep("Others",6),rep("Brain related",3),rep("Others",3),"Colon related", 
"Colon related", rep("Others",16),"Colon related","Others","Colon related",rep("Others",5))
gAnnotationData = data.frame(tissue_name, category)
gLogCpmData = a
pdf("Tissue_GTEx_heatmap_unscaled.pdf",width=8, height=8)
diag(gLogCpmData)=1    
testHeatmap3(gLogCpmData, gAnnotationData)
dev.off()




```

## QQ plot
https://uw-gac.github.io/topmed_workshop_2017/association-tests.html#association-testing-with-aggregate-units
```R
library(ggplot2)
qqPlot <- function(pval) {
    pval <- pval[!is.na(pval)]
    n <- length(pval)
    x <- 1:n
    dat <- data.frame(obs=sort(pval),
                      exp=x/n,
                      upper=qbeta(0.025, x, rev(x)),
                      lower=qbeta(0.975, x, rev(x)))
    
    ggplot(dat, aes(-log10(exp), -log10(obs))) +
        geom_line(aes(-log10(exp), -log10(upper)), color="gray") +
        geom_line(aes(-log10(exp), -log10(lower)), color="gray") +
        geom_point() +
        geom_abline(intercept=0, slope=1, color="red") +
        xlab(expression(paste(-log[10], "(expected P)"))) +
        ylab(expression(paste(-log[10], "(observed P)"))) +
        theme_bw()
}    

qqPlot(assoc$Wald.pval)
```

## QQplot multiple groups
I modified the codes from https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R

so that the points won't be along the y-axis when we observe many small p values 
```R
library(lattice)
qqunif.plot<-function(pvalues, 
	should.thin=T, thin.obs.places=2, thin.exp.places=2, 
	xlab=expression(paste("Expected (",-log[10], " p-value)")),
	ylab=expression(paste("Observed (",-log[10], " p-value)")), 
	draw.conf=TRUE, conf.points=1000, conf.col="lightgray", conf.alpha=.05,
	already.transformed=FALSE, pch=20, 
	aspect="iso", 
	prepanel=prepanel.qqunif,
	par.settings=list(superpose.symbol=list(pch=pch)), ...) {
	
	
	#error checking
	if (length(pvalues)==0) stop("pvalue vector is empty, can't draw plot")
	if(!(class(pvalues)=="numeric" || 
		(class(pvalues)=="list" && all(sapply(pvalues, class)=="numeric"))))
		stop("pvalue vector is not numeric, can't draw plot")
	if (any(is.na(unlist(pvalues)))) stop("pvalue vector contains NA values, can't draw plot")
	if (already.transformed==FALSE) {
		if (any(unlist(pvalues)==0)) stop("pvalue vector contains zeros, can't draw plot")
	} else {
		if (any(unlist(pvalues)<0)) stop("-log10 pvalue vector contains negative values, can't draw plot")
	}
	
	
	grp<-NULL
	n<-1
	exp.x<-c()
	if(is.list(pvalues)) {
		nn<-sapply(pvalues, length)
		rs<-cumsum(nn)
		re<-rs-nn+1
		n<-min(nn)
		if (!is.null(names(pvalues))) {
			grp=factor(rep(names(pvalues), nn), levels=names(pvalues))
			names(pvalues)<-NULL
		} else {
			grp=factor(rep(1:length(pvalues), nn))
		}
		pvo<-pvalues
		pvalues<-numeric(sum(nn))
		exp.x<-numeric(sum(nn))
		for(i in 1:length(pvo)) {
			if (!already.transformed) {
				pvalues[rs[i]:re[i]] <- -log10(pvo[[i]])
				exp.x[rs[i]:re[i]] <- -log10((rank(pvo[[i]], ties.method="first")-.5)/nn[i])
			} else {
				pvalues[rs[i]:re[i]] <- pvo[[i]]
				exp.x[rs[i]:re[i]] <- -log10((nn[i]+1-rank(pvo[[i]], ties.method="first")-.5)/(nn[i]+1))
			}
		}
	} else {
		n <- length(pvalues)+1
		if (!already.transformed) {
			exp.x <- -log10((rank(pvalues, ties.method="first")-.5)/n)
			pvalues <- -log10(pvalues)
		} else {
			exp.x <- -log10((n-rank(pvalues, ties.method="first")-.5)/n)
		}
	}


	#this is a helper function to draw the confidence interval
	panel.qqconf<-function(n, conf.points=1000, conf.col="gray", conf.alpha=.05, ...) {
		require(grid)
		conf.points = min(conf.points, n-1);
		mpts<-matrix(nrow=conf.points*2, ncol=2)
        	for(i in seq(from=1, to=conf.points)) {
            		mpts[i,1]<- -log10((i-.5)/n)
            		mpts[i,2]<- -log10(qbeta(1-conf.alpha/2, i, n-i))
            		mpts[conf.points*2+1-i,1]<- -log10((i-.5)/n)
            		mpts[conf.points*2+1-i,2]<- -log10(qbeta(conf.alpha/2, i, n-i))
        	}
        	grid.polygon(x=mpts[,1],y=mpts[,2], gp=gpar(fill=conf.col, lty=0), default.units="native")
    	}

	#reduce number of points to plot
	if (should.thin==T) {
		if (!is.null(grp)) {
			thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
				exp.x = round(exp.x, thin.exp.places),
				grp=grp))
			grp = thin$grp
		} else {
			thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
				exp.x = round(exp.x, thin.exp.places)))
		}
		pvalues <- thin$pvalues
		exp.x <- thin$exp.x
	}
	gc()
	
	prepanel.qqunif= function(x,y,...) {
		A = list()
		A$xlim = range(x)*1.02
		A$xlim[1]=0
		A$ylim = range(y)*1.02
		A$ylim[1]=0
		return(A)
	}

	#draw the plot
	xyplot(pvalues~exp.x, groups=grp, xlab=xlab, ylab=ylab, 
		#aspect=aspect,
		prepanel=prepanel, 
		scales=list(axs="i"), 
		pch=pch,
		panel = function(x, y, ...) {
			if (draw.conf) {
				panel.qqconf(n, conf.points=conf.points, 
					conf.col=conf.col, conf.alpha=conf.alpha)
			};
			panel.xyplot(x,y, ...);
			panel.abline(0,1);
		}, par.settings=par.settings, ...
	)
}

my.pvalue.list<-list("AA"=mydat[mydat$Population=="AA",1], "EA"=mydat[mydat$Population=="EA",1],
"AFA"=mydat[mydat$Population=="AFA",1],"CAU"=mydat[mydat$Population=="CAU",1],
"HIS"=mydat[mydat$Population=="HIS",1])

pdf(paste0("QQ_trait",traits[traitnum],"_test.pdf") ,width = 6, height = 6)
qqunif.plot(my.pvalue.list, auto.key=list(corner=c(.95,.05)))
dev.off()




```
