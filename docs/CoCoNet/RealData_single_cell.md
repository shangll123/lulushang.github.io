---
layout: default
title: RealData single cell
parent: CoCoNet
nav_order: 2
---


### Build cell type specific networks in GTEx single cell dataset. 

##### Set workpath

```
workpath = "~path/coconet_cell/panda"
```

##### Data downloaded from GTEx Portal:

```
wget https://storage.googleapis.com/gtex_additional_datasets/single_cell_data/GTEx_droncseq_hip_pcf.tar
```

##### load needed files
```
load("~path/outcome_cell_scale.RData")
gene_table = read.table("~path/genetable.txt", header=T)
motif = read.table("~path/pandas/PANDA_for_server/data/motif.txt") # same as in GTEx tissue PANDA input, downloaded from same website 
ppi = read.table("~path/pandas/PANDA_for_server/data/ppi.txt") # same as in GTEx tissue PANDA input, downloaded from same website 
expr = read.table("~path/single_cell_gtex/GTEx_droncseq_hip_pcf/GTEx_droncseq_hip_pcf.umi_counts.txt",header=T) # GTEx single cell expression
cell_anno = read.csv("~path/single_cell_gtex/GTEx_droncseq_hip_pcf/cell_annotation.csv",header=T) # cell type annotation for GTEx single cell expression

#> dim(expr)
#[1] 32111 14963
#> dim(cell_anno)
#[1] 14963     6
```

##### motif is a file of lots of TF-gene pairs, column1: TF names, column2: gene names 
##### the dimension is same as in the PANDA output file PANDA_tissues_regulatory.txt
##### PANDA_tissues_regulatory.txt contains all edge weights of TF-gene pairs in motif.txt

##### grch37 for converting HGNC to ENSG ID

```
library(biomaRt)
grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice")
listDatasets(grch37)
grch37 = useDataset("hsapiens_gene_ensembl",mart=grch37)
```
##### map HGNC ID to ENSG ID, since we are using ENSG IDs in gene effect size file and only keep genes that we already have annotation for TSS TES location
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
expr_gene = rownames(expr)[match(as.character(dictionary1$hgnc_symbol),rownames(expr) )]
expr_mat = expr[match(as.character(dictionary1$hgnc_symbol),rownames(expr) ),]
```

#### Prepare RPKM expression
```
gene_annotate = match(dictionary1$ensembl_gene_id,as.character(gene_table$ENSG_ID))
l = gene_table$distance[gene_annotate]
cS <- colSums(expr_mat) #Total mapped reads per sample
rpkm <- (10^6)*t(t(expr_mat/l)/cS)
#save(rpkm, file = "rpkm.RData")
```

##### double quantile Normalization on log10(rpkm+1)
```
log10_rpkm = log10(rpkm+1)			    		    
GTEx_expr_sc_log10_plus1_gene <- t(apply(log10_rpkm,2, function(x) qqnorm(x, plot=F)$x))		 
GTEx_expr_sc_log10_plus1_cell <- t(apply(GTEx_expr_sc_log10_plus1_gene,1, function(x) qqnorm(x, plot=F)$x))
GTEx_expr_sc_log10 = t(GTEx_expr_sc_log10_plus1_cell)
#save(GTEx_expr_sc_log10, file = "GTEx_expr_sc_log10.RData")				    
```
##### make cell type annotation. Get cell type index for each cell
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
library( matrixStats )
IQR_rpkmlog10_all = unlist(apply(GTEx_expr_sc_log10,1, function(x) IQR(x)))     
median_rpkmlog10_all = rowMedians(GTEx_expr_sc_log10) 
gene_rpkmlog10_specifity = matrix(0,19822,15)                      
for(i in 1:15){
	print(i)
	expr_tmp = GTEx_expr_sc_log10[, index_type[[i]]]
  	median_tmp = rowMedians(expr_tmp)
	gene_rpkmlog10_specifity[,i] = (median_tmp-median_rpkmlog10_all)/IQR_rpkmlog10_all
}                      
m = rowSums(gene_rpkmlog10_specifity)
```			    

##### retained genes with total specificity score in tissues greater than the median value across all genes. 

```				 
ind = which(m>-2.803) # use median
length(ind)
GTEx_expr_sc = GTEx_expr_sc_log10[ind,]
rownames(GTEx_expr_sc) = dictionary1$hgnc_symbol[ind]				 
dict9900 = dictionary1[ind,]
```

##### make new motif file. In the TF by gene matrix produced by PANDA, if std of any column is 0, it will cause NA in the normalization step, and algorithm won't work. in other words, we retained genes that are TF factors and have at least one connection to genes we already retained. 
```
genenames = rownames(GTEx_expr_sc)
motif_new = motif[which(as.character(motif$V1) %in% genenames),]			 
xx=intersect(as.character(motif_new$V2), dict9900$ensembl_gene_id)
motif_use = motif_new[which(as.character(motif_new$V2) %in% xx),]		 				 
dict8270 = dict9900[which(dict9900$ensembl_gene_id %in% xx),]		 
GTEx_expr_sc_8270 = GTEx_expr_sc[which(dict9900$ensembl_gene_id %in% xx),]	
```
##### need to remove gene with 0 std of i-th column in TF by gene matrix
```
i=6560					 				 
# update files		 
dict_update = dict8270[-6560,]	
GTEx_expr_sc_ensg = 	GTEx_expr_sc_8270[-6560,]
rownames(GTEx_expr_sc_ensg)=dict_update$ensembl_gene_id					 
motif_update = motif_use[-which(as.character(motif_use$V2)=="ENSG00000211746"),]		
```
##### make new ppi file

```			 
ppi[,1]=as.character(ppi[,1])
ppi[,2]=as.character(ppi[,2])				 
ind1 = which(ppi[,1] %in% dict_update$hgnc_symbol)
ind2 = which(ppi[,2] %in% dict_update$hgnc_symbol)		
ind = intersect(ind1, ind2)				 
ppi_new = ppi[ind,]				 
write.table(ppi_new, "ppi_new.txt",  col.names=F, quote=F,row.names=F)
```	

##### reorder cells, to prepare input expression file for PANDA

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
write.table(GTEx_expr_sc_reorder, "GTEx_expr_sc_reorder.txt",  col.names=F, quote=F)
write.table(motif_update, "motif_update.txt",  col.names=F, quote=F,row.names=F)
write.table(my_cell_type_reorder, "my_cell_type_reorder.txt",  col.names=F, quote=F,row.names=F)
```

##### run matlab
matlab codes can be found here: 
[here](https://drive.google.com/drive/folders/18sgoHHx_x03y6zNTSF69vz9hlmeStk-z?usp=sharing).

```
nohup matlab -nodisplay -nodesktop -nojvm -nosplash < RunPANDA.m &				 
```

##### then collect results
```	
% in matlab (faster than in R)
motif_file='motif_update.txt'; % motif prior
exp_file='GTEx_expr_sc_reorder.txt'; % expression data without headers
ppi_file='ppi_new.txt'; %  ppi prior
sample_file='my_cell_type_reorder.txt'; % information on sample order and tissues				 
fid=fopen(exp_file, 'r');
disp('fid=fopen(exp_file)')
headings=fgetl(fid);
NumConditions=length(regexp(headings, ' '));
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
##### go back to R

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
				 				 
```				 
	
