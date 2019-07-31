---
layout: default
title: Gene effect size
parent: CoCoNet
nav_order: 4
---

##### With SNP-level GWAS summary statistics and LD information from the reference panel, we obtained gene-level heritability estimates using [MQS](https://www.ncbi.nlm.nih.gov/pubmed/29515717). We scaled the gene-level heritability estimates by the number of SNPs in each gene. 

##### The sumstat.meta file could be downloaded from [here](https://drive.google.com/open?id=1GlTwFMafeB2k0bbvOJRD0ObGPcCSgmUj).
##### The genetable.txt could be downloaded from [here](https://drive.google.com/open?id=1XkyFp8_k1FLoYiaL_PYjYzusYoc8Lwz_).
```R
snp=read.csv("~/effectsize/data/sumstat.meta",header=T,sep=" ")
gene=read.csv("~/effectsize/data/genetable.txt",header=T,sep="\t")
#seperate snp and gene according to chromosome number
for (name in levels(as.factor(snp$CHR))){
  tmp=subset(snp,CHR==name)
  fn=paste('snp_chr_',gsub(' ','',name),".txt",sep='')
  write.csv(tmp,fn,row.names=FALSE)
}
for (name in levels(as.factor(gene$chr))){
  tmp=subset(gene,chr==name)
  fn=paste('gene_chr_',gsub(' ','',name),".txt",sep='')
  write.csv(tmp,fn,row.names=FALSE)
}
```


##### Get -/+1MB snps, snps in genes have no overlap
```R
chr_gene_index2=NULL
for(chr_id in 1:22) {

  gmap <- read.table(paste0("~/effectsize/Analysis/gene_snp_chr/gene_chr_", chr_id, ".txt"), sep = ",", header = T)
  smap <- read.table(paste0("~/effectsize/Analysis/gene_snp_chr/snp_chr_", chr_id, ".txt"), sep = ",", header = T)[, c(1, 2, 3)]
  num_gene <- dim(gmap)[1]
  snp_index <- NULL
  genmap <- NULL
  block_num <- 0
  count=0

  for(i in 1:num_gene) {

    snp_start <- gmap[i, 3] - 1000000
    snp_end <- gmap[i, 4] + 1000000
    dis_to_lastgene=abs(smap[, 3]-gmap[i-1, 3])
    dis_to_currentgene=abs(smap[, 3]-gmap[i, 3])
    dis_to_nextgene=abs(smap[, 3]-gmap[i+1, 3])
    
    temp_index <- which(smap[, 3] >= snp_start &  smap[, 3] <= snp_end & dis_to_currentgene<dis_to_lastgene &dis_to_currentgene<dis_to_nextgene)
    if(length(temp_index) < 20) {
      cat("There is no enough SNP in   ", i, "   gene", fill = T)
      next
    } else {
      count=count+1
      snp_index <- rbind(snp_index, temp_index[c(1, length(temp_index))])
      genmap <- rbind(genmap, gmap[i, ])
      block_num <- block_num + 1
      cat("There are",  length(temp_index), " SNPs in   ", i, "   gene", fill = T)
    }
    fn=paste0("~/effectsize/Analysis/snplist_gene/chr_", chr_id, "_gene_",count,".snpname", sep = "")
    chr_gene_index2=c(chr_gene_index2,i)
    #write.table(smap[temp_index,1], file=fn, quote = F, row.names = F, col.names = F)
	
  }
  colnames(snp_index) <- c("snp_start", "snp_end")
  write.table(snp_index, paste("~/effectsize/Analysis/OneMB_snp/chr_", chr_id, ".index", sep = ""), quote = F, row.names = F, col.names = T)
  write.table(genmap, paste("~/effectsize/Analysis/OneMB_snp/chr_", chr_id, ".genmap", sep = ""), quote = F, row.names = F, col.names = T)
  write.table(count,paste("~/effectsize/Analysis/OneMB_snp/chr_", chr_id, ".num", sep = ""), quote = F, row.names = F, col.names = F)
  }

```

##### get gene table
```R
chr_gene_index1=NULL
chr_gene_index2=NULL
for(i in 1:22){
chr_num=read.table(paste0("~/effectsize/Analysis/OneMB_snp/chr_", i, ".num"),header=F)$V1
chr_gene_index2=c(chr_gene_index2,1:chr_num)
chr_gene_index1=c(chr_gene_index1,rep(i,chr_num))
}

chrlist=list()
for(i in 1:22){
  chrlist[[i]]= read.table(paste0("~/effectsize/Analysis/OneMB_snp/chr_",i,".genmap"),header=T)
  chrgnum=read.table(paste0("~/effectsize/Analysis/OneMB_snp/chr_",i,".num"))$V1
  geneID=NULL
  for(j in 1:chrgnum){
    geneID[j]=strsplit(as.character(chrlist[[i]][,6]), "[.]")[[j]][1]
  }
  index_match=c(1:chrgnum)
  chrlist[[i]]=cbind(chrlist[[i]][,c(1,2,3,4)],geneID,index_match)
}
m=do.call(rbind,chrlist)
```


##### Build Genotype matrix:

```R
library(BEDMatrix)
chr_bed=NULL
for(i in 1:length(chr_gene_index1)){
 bed=BEDMatrix(paste0("~/effectsize/Analysis/plink_files/chr_",chr_gene_index1[i],"_gene_",chr_gene_index2[i],".bed"))
 chr_bed=cbind(chr_bed,bed[,])
}
dim(chr_bed)
```

##### Make plink usable
```
chmod u+x ~/effectsize/Analysis/plink2
```
##### use plink to extract SNPs corresponding to each gene
```
#!/bin/bash
#SBATCH --job-name=plink
#SBATCH --mem=20000
#SBATCH --array=1-22
#SBATCH --output=~/effectsize/Analysis/plink_files/out/plink%a.out
#SBATCH --error=~/effectsize/Analysis/plink_files/err/plink%a.err

bash

let k=0

for ((i=1; i<=1000; i++)); do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
Gfile=~/task1/refgeno/refEUR/ping_chr_${k}
gnum=`cat ~/effectsize/Analysis/OneMB_snp/chr_${i}.num`
for ((j=1; j<=${gnum}; j++)); do
snplist=~/effectsize/Analysis/snplist_gene/chr_${i}_gene_${j}.snpname
./plink2 --bfile ${Gfile} --extract ${snplist} --make-bed --out chr_${i}_gene_${j}
rm *log
rm *nosex
done
fi
done
```

##### Use gemma to get per-SNP heritability
```
# note:
# mkdir gemma_files/output/SCZ..., output files into seperate folders for 8 traits
```
gemma codes
```
#!/bin/bash
#SBATCH --job-name=goodluck_gemma
#SBATCH --mem=50000
#SBATCH --array=1-22
#SBATCH --output=~/effectsize/Analysis/gemma_files/out/gemma%a.out
#SBATCH --error=~/effectsize/Analysis/gemma_files/err/gemma%a.err

bash

GEMMA=~/effectsize/task1/SummaryData/gemma
INPATH1=~/effectsize/task1/sumstat
INPATH2=~/effectsize/Analysis/plink_files

let k=0
for ((i=1; i<=100; i++)); do
	let k=${k}+1
	if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
		for TRAIT in SCZ BIP BIPSCZ Alzheimer PBC UC CD IBD; do
				gnum=`cat ~/effectsize/Analysis/OneMB_snp/chr_${i}.num`
				BETAFILE=${INPATH1}/${TRAIT}.sumstats
			for ((j=1; j<=${gnum}; j++)); do
				BFILE=${INPATH2}/chr_${i}_gene_${j}
				${GEMMA} -beta ${BETAFILE} -bfile ${BFILE} -c pop.txt -vc 1 -o ${TRAIT}/chr_${i}_gene_${j}
done
done
fi
done
```

##### Extract per-SNP heritability

```R
Traits_name=c( "SCZ", "BIP", "BIPSCZ",  "Alzheimer","PBC", "CD", "UC", "IBD")

Sigma2_traits=matrix(0,49015,8)
count=0
for( k in Traits_name){
print(k)
count=count+1
sigma2=NULL
numchr=NULL
SIGMA2=NULL

for(i in 1:22){
			print(i)
  			numchr=read.table(paste0('~/effectsize/Analysis/OneMB_snp/chr_',i,'.num'))$V1
                       for(j in 1:numchr){
                         memo=paste0("~/effectsize/Analysis/gemma_files/output/",k,"/chr_",i,"_gene_",j)
                         gemma.res=read.delim(paste0(memo, ".log.txt", sep=""), head=FALSE)
                         sigma2=as.numeric(unlist(strsplit(as.character(gemma.res[17, ]), "  "))[2])
                         SIGMA2=c(SIGMA2,sigma2)
                       }
}
Sigma2_traits[,count]=SIGMA2
}
colnames(Sigma2_traits)=Traits_name
```

#want SNP number per gene and se(sigma2)
```R
Sigma2_traits=matrix(0,49015,8)
Sigma2_se_traits=matrix(0,49015,8)
PVE_Traits=matrix(0,49015,8)
PVE_se_Traits=matrix(0,49015,8)

Traits_name = c( "SCZ", "BIP", "BIPSCZ",  "Alzheimer","PBC", "CD", "UC", "IBD")
count=0
for( k in Traits_name){
print(k)
count=count+1
sigma2=NULL
sigma2_se=NULL
pve=NULL
pve_se=NULL
snpnum=NULL
numchr=NULL
SIGMA2=NULL
SIGMA2_SE=NULL
SNPNUM=NULL
PVE=NULL
PVE_SE=NULL
for(i in 1:22){
			print(i)
  			numchr=read.table(paste0('~/effectsize/Analysis/OneMB_snp/chr_',i,'.num'))$V1
                       for(j in 1:numchr){
                         memo=paste0("~/effectsize/Analysis/gemma_files/output/",k,"/chr_",i,"_gene_",j)
                         gemma.res=read.delim(paste0(memo, ".log.txt", sep=""), head=FALSE)
                         #sigma2=as.numeric(unlist(strsplit(as.character(gemma.res[17, ]), "  "))[2])
                         pve=as.numeric(unlist(strsplit(as.character(gemma.res[15, ]), "  "))[2])
                         pve_se=as.numeric(unlist(strsplit(as.character(gemma.res[16, ]), "  "))[2])
                         #sigma2_se=as.numeric(unlist(strsplit(as.character(gemma.res[18, ]), "  "))[2])
                         snpnum=read.table(paste0(memo, ".size.txt", sep=""), head=FALSE)$V1[1]
                         #SIGMA2=c(SIGMA2,sigma2)
                         #SIGMA2_SE=c(SIGMA2_SE,sigma2_se)
                         PVE=c(PVE,pve)
                         PVE_SE=c(PVE_SE,pve_se)
                         SNPNUM=c(SNPNUM,snpnum)
                       }
}
PVE_Traits[,count]=PVE
PVE_se_Traits[,count]=PVE_SE
}
colnames(PVE_Traits)=Traits_name
colnames(PVE_se_Traits)=Traits_name

rownames(PVE_Traits)=m$geneID
rownames(PVE_se_Traits)=m$geneID
```
