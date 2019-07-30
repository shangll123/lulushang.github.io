---
layout: default
title: RealData single cell
parent: CoCoNet
nav_order: 2
---



<img align="left" src="/images/coconuts.png" alt="drawing" width="50"/>

#### Set workpath

```
workpath = "/home/lulushang/pairwise_project/coconet_cell/panda"
```

#### Data downloaded from GTEx Portal:

```
wget https://storage.googleapis.com/gtex_additional_datasets/single_cell_data/GTEx_droncseq_hip_pcf.tar
```

#### load needed files
```
load("/home/lulushang/pairwise_project/GTEx/data/Sigma_pop_43traits.RData")
gene_table = read.table("/home/lulushang/pairwise_project/GTEx/data/genetable.txt", header=T)
motif = read.table("/home/lulushang/pairwise_project/pandas/PANDA_for_server/data/motif.txt") # same as in GTEx tissue PANDA input, downloaded from same website 
ppi = read.table("/home/lulushang/pairwise_project/pandas/PANDA_for_server/data/ppi.txt") # same as in GTEx tissue PANDA input, downloaded from same website 
expr = read.table("/home/lulushang/pairwise_project/single_cell_gtex/GTEx_droncseq_hip_pcf/GTEx_droncseq_hip_pcf.umi_counts.txt",header=T) # GTEx single cell expression
cell_anno = read.csv("/home/lulushang/pairwise_project/single_cell_gtex/GTEx_droncseq_hip_pcf/cell_annotation.csv",header=T) # cell type annotation for GTEx single cell expression


#> dim(expr)
#[1] 32111 14963
#> dim(cell_anno)
#[1] 14963     6

# motif is a file of lots of TF-gene pairs, column1: TF names, column2: gene names
# the dimension is same as in the PANDA output file PANDA_tissues_regulatory.txt
# PANDA_tissues_regulatory.txt contains all edge weights of TF-gene pairs in motif.txt
```


#### grch37 for converting HGNC to ENSG ID

```
library(biomaRt)
grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice")
listDatasets(grch37)
grch37 = useDataset("hsapiens_gene_ensembl",mart=grch37)
```
#### map HGNC ID to ENSG ID, since we are using ENSG IDs in gene effect size file and only keep genes that we already have annotation for TSS TES location
```
dict = rownames(expr)
HGNC_to_ENSE=getBM(attributes=c('hgnc_symbol','ensembl_gene_id'), 
       filters = 'hgnc_symbol', 
       values = dict, 
       mart = grch37)  # convert gene from  hgnc_symbol to ensembl_gene_id, make dictionary
       
ENSG_gene_table = as.character(gene_table$ENSG_ID)
index = which(HGNC_to_ENSE$ensembl_gene_id %in% ENSG_gene_table ) 
dictionary = HGNC_to_ENSE[index,]
dictionary1 = dictionary[-13896,] # the row 13896 is empty

dim(dictionary1)
head(dictionary1)

#> dim(dictionary1)
#[1] 19822     2
#> head(dictionary1)
#   hgnc_symbol ensembl_gene_id
#1         A1BG ENSG00000121410
#3         A1CF ENSG00000148584
#4          A2M ENSG00000175899
#5      A2M-AS1 ENSG00000245105
#9      A3GALT2 ENSG00000184389
#10      A4GALT ENSG00000128274


expr_gene = rownames(expr)[match(as.character(dictionary1$hgnc_symbol),rownames(expr) )]
expr_mat = expr[match(as.character(dictionary1$hgnc_symbol),rownames(expr) ),]

#> sum(rownames(expr_mat)==dictionary1$hgnc_symbol)
#[1] 19822
#> dim(expr_mat)
#[1] 19822 14963
# the gene name by order in dictionary1 and expr_mat are matched 

#save(dictionary1, file = "dictionary1.RData")
#save(expr_mat, file = "expr_mat.RData")
```
#### Prepare RPKM expression
```
gene_annotate = match(dictionary1$ensembl_gene_id,as.character(gene_table$ENSG_ID))
l = gene_table$distance[gene_annotate]
cS <- colSums(expr_mat) #Total mapped reads per sample
rpkm <- (10^6)*t(t(expr_mat/l)/cS)
#save(rpkm, file = "rpkm.RData")
```
#### double quantile Normalization on log10(rpkm+1)
```
log10_rpkm = log10(rpkm+1)			    
save(log10_rpkm, file = "log10_rpkm.RData")			    
			    
GTEx_expr_sc_log10_plus1_gene <- t(apply(log10_rpkm,2, function(x) qqnorm(x, plot=F)$x))		 
GTEx_expr_sc_log10_plus1_cell <- t(apply(GTEx_expr_sc_log10_plus1_gene,1, function(x) qqnorm(x, plot=F)$x))
GTEx_expr_sc_log10 = t(GTEx_expr_sc_log10_plus1_cell)
#save(GTEx_expr_sc_log10, file = "GTEx_expr_sc_log10.RData")				    
```


#### make cell type annotation. Get cell type index for each cell
```
sum(colnames(expr_norm)==as.character(cell_anno$Cell_ID)) 
colnameexpr = gsub(".", '-', as.character(colnames(expr_norm)), fixed = T)
type_in_order = c()
for(i in 1:length(colnames(expr_norm))){
	type_in_order[i] = which(as.character(cell_anno$Cell_ID) %in% colnameexpr[i])
}
cell_type = as.character(cell_anno$Cluster_Name)[type_in_order]
celluniquetype = unique(cell_type)
index_type  = list()
count = 0
for(i in celluniquetype){
	count = count + 1
	index_type[[count]] = which(cell_type == i)
}
```	 

#### calculate gene specificity by cell type					 
```					 
IQR_rpkmlog10_all = unlist(apply(GTEx_expr_sc_log10,1, function(x) IQR(x)))
        
library( matrixStats )
median_rpkmlog10_all = rowMedians(GTEx_expr_sc_log10) 
gene_rpkmlog10_specifity = matrix(0,19822,15)                      
for(i in 1:15){
	print(i)
	expr_tmp = GTEx_expr_sc_log10[, index_type[[i]]]
  	median_tmp = rowMedians(expr_tmp)
	gene_rpkmlog10_specifity[,i] = (median_tmp-median_rpkmlog10_all)/IQR_rpkmlog10_all
}
                       
m = rowSums(gene_rpkmlog10_specifity)
			    
#> summary(m)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#-20.587  -2.953  -2.803  -1.726  -2.663 794.171      14 

# retained genes with total specificity score in tissues greater than the median value across all genes. 

				 
ind = which(m>-2.803) # use median
length(ind)
GTEx_expr_sc = GTEx_expr_sc_log10[ind,]
rownames(GTEx_expr_sc) = dictionary1$hgnc_symbol[ind]	
				 
dict9900 = dictionary1[ind,]
#save(dict9900,file = "dict9900.RData")			 
#save(GTEx_expr_sc, file = "GTEx_expr_sc.RData")	
```

#### make new motif file. In the TF by gene matrix produced by PANDA, if std of any column is 0, it will cause NA in the normalization step, and algorithm won't work. in other words, we retained genes that are TF factors and have at least one connection to genes we already retained. 
```
genenames = rownames(GTEx_expr_sc)
motif_new = motif[which(as.character(motif$V1) %in% genenames),]
#> dim(motif_new)
#[1] 7832937       3
#> dim(motif)
#[1] 19476492        3				 
```
#### check if genes in dict9900 could be found in motif_new
```
#> ind = which(as.character(motif_new$V1) %in% dict9900$hgnc_symbol)
#> length(ind)
#[1] 7832937
#> TF=unique(as.character(motif_new$V1))
#> length(TF)
#[1] 259
#> length(intersect(TF, dict9900$hgnc_symbol))
#[1] 259
#> length(intersect(as.character(motif_new$V2), dict9900$ensembl_gene_id))
#[1] 8270
xx=intersect(as.character(motif_new$V2), dict9900$ensembl_gene_id)
motif_use = motif_new[which(as.character(motif_new$V2) %in% xx),]		 				 
#> dim(motif_use)
#[1] 2141930       3
				 
dict8270 = dict9900[which(dict9900$ensembl_gene_id %in% xx),]
# save(	dict8270, file = "dict8270.RData")			 
GTEx_expr_sc_8270 = GTEx_expr_sc[which(dict9900$ensembl_gene_id %in% xx),]	
```
#### need to remove gene with 0 std of i-th column in TF by gene matrix
```
i=6560	
				 
#> dict8270[6560,]
#      hgnc_symbol ensembl_gene_id
#21956      TRBV19 ENSG00000211746	
				 
# update files
			 
dict_update = dict8270[-6560,]	
GTEx_expr_sc_ensg = 	GTEx_expr_sc_8270[-6560,]
rownames(GTEx_expr_sc_ensg)=dict_update$ensembl_gene_id					 
motif_update = motif_use[-which(as.character(motif_use$V2)=="ENSG00000211746"),]

#> dim(motif_update)
# [1] 2141671       3	
#> dict8269 = dict8270[-6560,]
#> save(dict8269,file = "dict8269.RData")		

#-------------------------
# make new ppi file
#-------------------------				 

ppi[,1]=as.character(ppi[,1])
ppi[,2]=as.character(ppi[,2])				 
ind1 = which(ppi[,1] %in% dict_update$hgnc_symbol)
ind2 = which(ppi[,2] %in% dict_update$hgnc_symbol)		
ind = intersect(ind1, ind2)
#> length(ind)
#[1] 13883				 
ppi_new = ppi[ind,]				 
write.table(ppi_new, "ppi_new.txt",  col.names=F, quote=F,row.names=F)
```		 	 
#### reorder cells, to prepare input expression file for PANDA
```
tissue_names = c( "ASC1"     ,    "ASC2"   ,      "END"       ,   "GABA1"    ,    "GABA2"       ,
"MG"     ,      "NSC"        ,  "ODC1"       ,  "OPC"    ,      "Unclassified",
"exCA1"     ,   "exCA3"  ,      "exDG"     ,    "exPFC1"    ,   "exPFC2" )

expr_list = list()
count = 0
GTEx_expr_sc_reorder = matrix(0,dim(GTEx_expr_sc_ensg)[1],1)
cell_type_reorder = c()
for(i in tissue_names){
	print(count)
	count = count+1
	col_index = which(cell_type %in% i)
	expr_list[[count]] = GTEx_expr_sc_ensg[,col_index]
	GTEx_expr_sc_reorder = cbind(GTEx_expr_sc_reorder,expr_list[[count]] )
	cell_type_reorder = c(cell_type_reorder,rep(i,dim(expr_list[[count]])[2]))
}
GTEx_expr_sc_reorder = GTEx_expr_sc_reorder[,-1]
my_cell_type_reorder = cbind(colnames(GTEx_expr_sc_reorder), cell_type_reorder)
				 
#save(GTEx_expr_sc_reorder, file = "GTEx_expr_sc_reorder.RData")
#save(motif_new, file = "motif_new.RData")
#save(my_cell_type_reorder, file = "my_cell_type_reorder.RData")	
			 
write.table(GTEx_expr_sc_reorder, "GTEx_expr_sc_reorder.txt",  col.names=F, quote=F)
write.table(motif_update, "motif_update.txt",  col.names=F, quote=F,row.names=F)
write.table(my_cell_type_reorder, "my_cell_type_reorder.txt",  col.names=F, quote=F,row.names=F)
```

#### run matlab
```
nohup matlab -nodisplay -nodesktop -nojvm -nosplash < RunPANDA.m &				 
```

		 
 
#### then collect results
```		 
motif_file='motif_update.txt'; % motif prior
exp_file='GTEx_expr_sc_reorder.txt'; % expression data without headers
ppi_file='ppi_new.txt'; %  ppi prior
sample_file='my_cell_type_reorder.txt'; % information on sample order and tissues
% 
% 				 
fid=fopen(exp_file, 'r');
disp('fid=fopen(exp_file)')
headings=fgetl(fid);
NumConditions=length(regexp(headings, ' '));
% 
% 
disp('NumConditions')
frewind(fid);
Exp=textscan(fid, ['%s', repmat('%f', 1, NumConditions)], 'delimiter', ' ', 'CommentStyle', '#');
fclose(fid);
GeneNames=Exp{1};
disp('GeneNames')
NumGenes=length(GeneNames);
disp('NumGenes')
Exp=cat(2, Exp{2:end});
disp('Exp')
% 
% 			 
PANDA_InputEdges = 'PANDA_InputEdges.pairs'; 
PANDA_weights = 'PANDA_tissues_regulatory.txt'				 
[TF, gene, weight]=textread(PANDA_InputEdges, '%s%s%f');
disp('PANDA_InputEdges read')
[weight1, weight2, weight3,weight4, weight5, weight6,weight7, weight8, weight9,weight10, weight11, weight12,weight13, weight14, weight15]=textread(PANDA_weights, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f ');				 
disp('PANDA_weights read')
TFNames=unique(TF);
disp('TFNames')
NumTFs=length(TFNames);
[~,i]=ismember(TF, TFNames);
[~,j]=ismember(gene, GeneNames);
%				 
RegNet=zeros(NumTFs, NumGenes);
RegNet(sub2ind([NumTFs, NumGenes], i, j))=weight1;
save('RegNet1.mat','RegNet')				 
% do above for 1-15 tissues	
				 
save('TF_gene_names.mat','TFNames','GeneNames')
```


#### go back to R				 
```			 
library(R.matlab)
TF_gene_names = readMat("TF_gene_names.mat")					 
i=1
for(i in 1:15){	
print(i)	
Test_data = readMat(paste0("RegNet",i,".mat"))				 
TF_gene_mat = Test_data$RegNet			 
colnames(TF_gene_mat) = unlist(TF_gene_names$GeneNames)	
rownames(TF_gene_mat) = unlist(TF_gene_names$TFNames)	# gene symbol				 
rownames(TF_gene_mat) = dict8269$ensembl_gene_id[which(dict8269$hgnc_symbol %in% rownames(TF_gene_mat))]# gene ENSG id				 
mat1 = matrix(0,dim(TF_gene_mat)[2],dim(TF_gene_mat)[2])
mat2 = matrix(0,dim(TF_gene_mat)[2],dim(TF_gene_mat)[2])				 
TF_gene_mat = (TF_gene_mat>0)		 
mat1[which(colnames(TF_gene_mat) %in% rownames(TF_gene_mat)),]=TF_gene_mat		
mat2[,which(colnames(TF_gene_mat) %in% rownames(TF_gene_mat))]=t(TF_gene_mat)	
mat = mat1+mat2	
mat[mat>0]=1				 
rownames(mat) = colnames(TF_gene_mat)				 
colnames(mat) = colnames(TF_gene_mat)					 
save(mat, file = paste0("Coconet_cell_",i,".RData"))
}				 
				 
#test=TF_gene_mat[,which(colnames(TF_gene_mat) %in% rownames(TF_gene_mat))]			 
#> test[1:5,1:5]
#                ENSG00000106546 ENSG00000116017 ENSG00000150347 ENSG00000143437
#ENSG00000106546        -0.38413         5.11820       -0.303530        -0.49125
#ENSG00000116017         1.67850         1.92720        2.684100         1.32500
#ENSG00000150347        -0.15517        -0.71310        0.306340         0.53523
#ENSG00000143437        -0.40894        -0.31624        0.047128        -0.25759				 
	
#---------------------------------------
				 
load("/home/lulushang/pairwise_project/Sigma_pop_43traits.RData")				 
index = which(as.character(gene_table$ENSG_ID) %in% dict8269$ensembl_gene_id)
sigma_cell = Sigma_pop_43traits[index,]
#save(sigma_cell, file = "sigma_cell.RData")
```				 
				 
				 



 
