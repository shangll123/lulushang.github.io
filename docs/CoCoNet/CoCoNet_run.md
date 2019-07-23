---
layout: default
title: CoCoNet
parent: CoCoNet
nav_order: 1
---

<!--- [_config.yml]({{ site.baseurl }}/images/config.png)--->

<!--- <img src="https://latex.codecogs.com/svg.latex?\Large&space;x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" title="\Large x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" />--->

<!--- <img src="https://latex.codecogs.com/svg.latex?\sum&space;\bigcup_{1}^{n}\overleftarrow{abc}" title="\sum \bigcup_{1}^{n}\overleftarrow{abc}" /> --->



#### Our Paper on bioRxiv: [Leveraging Gene Co-expression Patterns to Infer Trait-Relevant Tissues in Genome-wide Association Studies](https://www.biorxiv.org/content/biorxiv/early/2019/07/17/705129.full.pdf)

CoCoNet incorporates tissue-specific gene co-expression networks constructed from either bulk or single cell RNA sequencing (RNAseq) studies with GWAS data for trait-tissue inference. In particular, CoCoNet relies on a covariance regression network model to express gene-level effect sizes for the given GWAS trait as a function of the tissue-specific co-expression adjacency matrix. With a composite likelihood-based inference algorithm, CoCoNet is scalable to tens of thousands of genes. 
 

## Schematic 

![CoCoNet\ Schematic Overview](https://raw.githubusercontent.com/shangll123/shangll123.github.io/master/images/CoCoNet_Figure/Figure1.tiff)

**CoCoNet** is an efficient method to facilitate the identification of trait-relevant tissues or cell types. We apply CoCoNet for an in-depth analysis of four neurological disorders and four autoimmune diseases, where we integrate the corresponding GWASs with bulk RNAseq data from 38 tissues and single cell RNAseq data from 10 cell types. 

## Install the Package
```R
# Install devtools if necessary
install.packages('devtools')

# Install CoCoNet
devtools::install_github('shangll123/CoCoNet')

# Load the package
library(CoCoNet)

```

## Example: GTEx tissues in GWAS trait BIPSCZ

Load the `CoCoNet` package and data, which can be downloaded from this google drive [here](https://drive.google.com/open?id=1XkyFp8_k1FLoYiaL_PYjYzusYoc8Lwz_).
```R
load("tissue_net.RData")
load("tissue_name.RData")
load("outcome_tissue_scale.RData")
```

In total we have 38 tissues, the network are ordered by the tissue names
```R
> tissue_name
 [1] "Adipose_subcutaneous"      "Adipose_visceral"         
 [3] "Adrenal_gland"             "Artery_aorta"             
 [5] "Artery_coronary"           "Artery_tibial"            
 [7] "Brain_other"               "Brain_cerebellum"         
 [9] "Brain_basal_ganglia"       "Breast"                   
[11] "Lymphoblastoid_cell_line"  "Fibroblast_cell_line"     
[13] "Colon_sigmoid"             "Colon_transverse"         
[15] "Gastroesophageal_junction" "Esophagus_mucosa"         
[17] "Esophagus_muscularis"      "Heart_atrial_appendage"   
[19] "Heart_left_ventricle"      "Kidney_cortex"            
[21] "Liver"                     "Lung"                     
[23] "Minor_salivary_gland"      "Skeletal_muscle"          
[25] "Tibial_nerve"              "Ovary"                    
[27] "Pancreas"                  "Pituitary"                
[29] "Prostate"                  "Skin"                     
[31] "Intestine_terminal_ileum"  "Spleen"                   
[33] "Stomach"                   "Testis"                   
[35] "Thyroid"                   "Uterus"                   
[37] "Vagina"                    "Whole_blood"  
```
The first tissue is "Adipose_subcutaneous", which looks like this:
```
> tissue_net[[1]][1:4,1:4]
                ENSG00000106546 ENSG00000160224 ENSG00000156150 ENSG00000052850
ENSG00000106546               0               0               0               0
ENSG00000160224               0               0               0               0
ENSG00000156150               0               0               0               0
ENSG00000052850               0               0               0               0
```


