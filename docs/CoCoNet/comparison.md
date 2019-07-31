---
layout: default
title: Comparison
parent: CoCoNet
nav_order: 5
---


##### We partially validated the identified trait-relevant tissue/cell types for the GWAS diseases by searching PubMed, by following [this paper](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007186). We reasoned that, if the tissue or cell type is relevant to the disease of interest, then there would be previous publications studying the disease in the particular tissue or cell type. Therefore, by counting the number of previous publications using the key word pairs of trait and tissue or the key word pairs of trait and cell type, we would have a reasonable quantitative estimate on the relevance between the trait and the corresponding tissue/cell type. 

We used an R package [RISmed](https://cran.r-project.org/web/packages/RISmed/index.html) to efficiently count the number of publications in PubMed that contain the names of both the trait and the tissue/cell type either in the abstract or in the title.

```R

library(RISmed)

#-------------------------------
#  GTEx Tissue
#-------------------------------

tissues_name = c(
"adipose subcutaneous",      "adipose visceral",         
"adrenal gland"  ,           "artery aorta"  ,           
"artery coronary",           "artery tibial"     ,       
"brain other",               "brain cerebellum"  ,       
"brain basal ganglia"    ,   "breast"    ,               
"lymphoblastoid cell line" , "fibroblast cell line"   ,  
"colon sigmoid"      ,       "colon transverse"  ,       
"gastroesophageal junction", "esophagus mucosa"    ,     
"esophagus muscularis",      "heart atrial appendage"   ,
"heart left ventricle" ,     "kidney cortex"     ,       
"liver"          ,           "lung"       ,              
"minor salivary gland"  ,    "skeletal muscle"     ,     
"tibial nerve" ,             "ovary"   ,                 
"pancreas"  ,                "pituitary"       ,         
"prostate"  ,                "skin"    ,                 
"intestine terminal ileum" , "spleen"       ,            
"stomach"  ,                 "testis"   ,                
"thyroid"  ,                 "uterus"      ,             
"vagina"   ,                 "whole blood"   
)

trait_name = c(
"Schizophrenia",
"Bipolar disorder",
"Bipolar disorder & Schizophrenia",
"Alzheimer's disease",
"Primary Biliary Cholangitis",
"Crohn's disease",
"Ulcerative colitis",
"Inflammatory bowel disease"
)

tissues = c(
"(adipose subcutaneous[Title/Abstract] OR subcutaneous adipose[Title/Abstract])",  
"(adipose visceral[Title/Abstract] OR visceral adipose[Title/Abstract])",      
"(adrenal gland[Title/Abstract])"  ,          
"artery aorta[Title/Abstract]"  ,           
"artery coronary[Title/Abstract]",           
"artery tibial[Title/Abstract]",       
"(substantia nigra[Title/Abstract] OR hypothalamus[Title/Abstract] OR hippocampus[Title/Abstract] OR frontal lobe[Title/Abstract] OR cerebral cortex[Title/Abstract] OR amygdala[Title/Abstract])",               		 
"(cerebellum[Title/Abstract])",       
"(nucleus accumbens[Title/Abstract] OR caudate putamen[Title/Abstract] OR caudate nucleus[Title/Abstract])"    ,   
"breast[Title/Abstract]", 
"(lymphoblastoid cell line[Title/Abstract])" , 
"(fibroblast cell line[Title/Abstract])"   ,  
"(colon sigmoid[Title/Abstract] OR large intestine[Title/Abstract])"      ,
"(colon transverse[Title/Abstract] OR large intestine[Title/Abstract])"  ,       
"(gastroesophageal junction[Title/Abstract])",
"(esophagus mucosa[Title/Abstract])"    ,     
"(esophagus muscularis[Title/Abstract])",     
"(heart atrial appendage[Title/Abstract])"   ,
"(heart left ventricle[Title/Abstract])" ,
"(kidney cortex[Title/Abstract])"     ,       
"liver[Title/Abstract]",       
"lung[Title/Abstract]",              
"(minor salivary gland[Title/Abstract])"  ,   
"(skeletal muscle[Title/Abstract])"     ,     
"(tibial nerve[Title/Abstract])" ,            
"ovary[Title/Abstract]"   ,                 
"pancreas[Title/Abstract]"  ,                
"pituitary[Title/Abstract]"       ,         
"prostate[Title/Abstract]"  ,               
"skin[Title/Abstract]"    ,                 
"(intestine terminal ileum[Title/Abstract] OR terminal ileum[Title/Abstract])" , 
"spleen[Title/Abstract]"       ,            
"stomach[Title/Abstract]"  ,                 
"testis[Title/Abstract]"   ,                
"thyroid[Title/Abstract]"  ,                 
"uterus[Title/Abstract]"      ,             
"vagina[Title/Abstract]"   ,                 
"whole blood[Title/Abstract]"   
)


traits = c(
"Schizophrenia[Title/Abstract]",
"Bipolar disorder[Title/Abstract]",
"(Bipolar disorder Schizophrenia[Title/Abstract] OR Bipolar disorder[Title/Abstract] OR Schizophrenia[Title/Abstract])",
"(Alzheimer's disease[Title/Abstract] OR Alzheimer[Title/Abstract])",
"(Primary Biliary Cholangitis[Title/Abstract] OR PBC[Title/Abstract])",
"(Crohn's disease[Title/Abstract] OR Crohn[Title/Abstract])",
"(Ulcerative colitis[Title/Abstract] OR Ulcerative[Title/Abstract])",
"(Inflammatory bowel[Title/Abstract] OR IBD[Title/Abstract])"
)

results_all = matrix(0,8,38)
for(i in 1:8){
	print(i)
	for(j in 1:38){
		print(j)
		topic <- paste0(traits[i]," AND ",tissues[j])
		r <- QueryCount(EUtilsSummary(topic, db= "pubmed"))
		results_all[i,j] = r
		Sys.sleep(0.1)
}
}
colnames(results_all) = tissues_name
rownames(results_all) = trait_name

write.csv(results_all, "results_tissue.csv",quote=F)
results_all_percentage = t(apply(results_all,1,function(x) x/sum(x)))
results_all_percentage_use = round(results_all_percentage,4)*100
write.csv(results_all_percentage_use, "pubmed_result_tissue.csv",quote=F)

#-------------------------------
# GTEx Cell
#-------------------------------

cellnames = c("ASC","END","GABA","MG","NSC","ODC","OPC","exCA","exDG","exPFC")
traits = c(
"Schizophrenia[Title/Abstract]",
"Bipolar disorder[Title/Abstract]",
"(Bipolar disorder Schizophrenia[Title/Abstract] OR Bipolar disorder[Title/Abstract] OR Schizophrenia[Title/Abstract])",
"(Alzheimer's disease[Title/Abstract] OR Alzheimer[Title/Abstract])",
"(Primary Biliary Cholangitis[Title/Abstract] OR PBC[Title/Abstract])",
"(Crohn's disease[Title/Abstract] OR Crohn[Title/Abstract])",
"(Ulcerative colitis[Title/Abstract] OR Ulcerative[Title/Abstract])",
"(Inflammatory bowel[Title/Abstract] OR IBD[Title/Abstract])"
)

celltype = c(
"astrocytes[Title/Abstract]",
"endothelial[Title/Abstract]",
"GABAergic[Title/Abstract]",
"microglia[Title/Abstract]",
"neuronal stem cells[Title/Abstract]",
"oligodendrocytes[Title/Abstract]",
"oligodendrocyte precursor cells[Title/Abstract]",
"pyramidal neurons[Title/Abstract]",
"granule neurons[Title/Abstract]",
"glutamatergic neurons[Title/Abstract]"
)

results_all_cell = matrix(0,8,10)
for(i in 1:8){
	print(i)
	for(j in 1:10){
		topic <- paste0(traits[i]," AND ",celltype[j])
		r <- QueryCount(EUtilsSummary(topic, db= "pubmed"))
		results_all_cell[i,j] = r
		Sys.sleep(0.1)
}
}
colnames(results_all_cell) = cellnames
rownames(results_all_cell) = trait_name
write.csv(results_all_cell, "results_celltype.csv",quote=F)

results_all_percentage_cell = t(apply(results_all_cell,1,function(x) x/sum(x)))
results_all_percentage_use_cell = round(results_all_percentage_cell,4)*100
write.csv(results_all_percentage_use_cell, "pubmed_result_celltype.csv",quote=F)




```
