---
layout: default
title: Network - Tissue
parent: CoCoNet
nav_order: 3
---


For the tissue-specific network construction, we used the result from this paper:
[Understanding Tissue-Specific Gene Regulation](https://www.sciencedirect.com/science/article/pii/S2211124717314183?via%3Dihub)


This paper reconstructed networks from panda, built Gene regulatory networks for 38 human tissues. All needed data can be found [here](https://zenodo.org/record/838734#.XB0xoy3MwWo).


##### First check how many genes in GTEx_PANDA_tissues.RData overlap with our gene level effect sizes extracted from summary statistics.
```
load("GTEx_PANDA_tissues.RData") 
# annotate genes with grch37 
library(biomaRt)
grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice")
listDatasets(grch37)
grch37 = useDataset("hsapiens_gene_ensembl",mart=grch37)
# extract genes overlapped with 49015 genes in summary statistics
load("~/path/Sigma_pop_traits.RData")
gene_table = read.table("~/path/genetable.txt", header=T)
TF = as.character(edges$TF)
dict = unique(TF) # make dictionary
HGNC_to_ENSE=getBM(attributes=c('hgnc_symbol','ensembl_gene_id'), 
       filters = 'hgnc_symbol', 
       values = dict, 
       mart = grch37)
ENSG_gene_table = as.character(gene_table$ENSG_ID)
index = which(HGNC_to_ENSE$ensembl_gene_id %in% ENSG_gene_table )
dict = HGNC_to_ENSE[index,]
```
##### get edges with both genes exist in ENSG_gene_table
```
# filter genes in first column
index_col_1 = which(as.character(edges$TF) %in% as.character(dict$hgnc_symbol) )
# filter genes in second column
ENSG_gene_table = as.character(gene_table$ENSG_ID)
index_col_2 = which(as.character(edges$Gene) %in% ENSG_gene_table )
index_col_edge = intersect(index_col_1, index_col_2)
#> length(unique(as.character(edges$Gene)[index_col_edge]))
#[1] 25991
```

##### filter genes by specificity in tissues
```
net_specific = rowSums(netTS)
index_specific = which(net_specific>0)
index_specific_col = intersect(index_specific, index_col_edge)
exp_specific = rowSums(expTS)
gene_expr_ge = as.character(genes$Name[which(exp_specific>1)])
index_gene_expr_ge = which(as.character(edges$Gene) %in% gene_expr_ge)
index_specific = which(net_specific>0)
index_specific_col = intersect(index_specific, index_col_edge)
ind = intersect(index_gene_expr_ge, index_specific_col)
use_gene = unique(as.character(edges$Gene)[ind])
use_edges = edges[ind,]
use_net = net[ind,]
all_gene = unique(c(as.character(dict$ensembl_gene_id),as.character(use_edges$Gene)))

#> length(all_gene)
#[1] 5359

save(all_gene, file = "all_gene.RData")
save(use_net,file = "use_net.RData")
save(use_edges,file = "use_edges.RData")
```

##### get gene effect sizes for 5359 genes
```
tissue_name = colnames(use_net)
save(tissue_name, file = "tissue_name.RData")
gene_tsn = all_gene
save(gene_tsn, file = "gene_tsn.RData")
ENSG_gene_table = as.character(gene_table$ENSG_ID)
indexk = NULL
for(k in 1:length(all_gene)){
	indexk[k] = which(ENSG_gene_table %in% all_gene[k])
}

sigma_sona = Sigma_pop_traits[indexk,]
rownames(sigma_sona) = gene_tsn
save(sigma_sona, file = "outcome_tissue.RData")

```
##### Build network using edge information

 
```

# this part of R code can be slow, better to refer to the single cell process procedures and use matlab

name1 = as.character(use_edges$TF)
name2 = as.character(use_edges$Gene)
index1 = NULL
index2 = NULL
weightindex = NULL
for(i in 1:length(name1)){
if(length(which(dict$hgnc_symbol %in% name1[i]))>0){
name1_ensg = as.character(dict$ensembl_gene_id)[which(dict$hgnc_symbol %in% name1[i])]
index1[i] = which(all_gene %in% name1_ensg)
index2[i] = which(all_gene %in% name2[i])
}
}
index_dat = data.frame(index1, index2)
save(index_dat, file = "index_dat.RData")

for(j in 1:38){
print(j)
mat = matrix(0, length(all_gene), length(all_gene))
for(xx in 1:length(index1)){
	mat[index1[xx],index2[xx]]=use_net[xx,j]
}
save(mat, file = paste0("mat_weighted_tissue_",j,".RData"))

}
```

