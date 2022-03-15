---
layout: default
title: Simulation
nav_order: 6
has_children: false
parent: SpatialPCA
permalink: /docs/Projects/SpatialPCA/Simulation
---

##### 
```R
simu = function(
  location,
  label,
  init_params,
  scenario,
  J,
  batch_facLoc,
  de_prop,
  de_facLoc,
  de_facScale,
  sim_seed,
  debug = FALSE
  )
{ 

  N <- nrow(location)
  set.seed(sim_seed)
  # 1.simulate count data
  noBatch <- ifelse(batch_facLoc == 0, TRUE, FALSE)
  params <- setParams(
    init_params,
    batchCells = rep(3 * N, 1),
    batch.rmEffect = noBatch,
    batch.facLoc = batch_facLoc,
    nGenes = J,
    group.prob = c(0.25, 0.25, 0.25,0.25),
    out.prob = 0,
    de.prob = de_prop,
    de.facLoc = de_facLoc,
    de.facScale = de_facScale,
    seed = sim_seed)

  sim_groups <- splatSimulate(
    params = params,
    method = "groups",
    verbose = FALSE)

  # remove cells having no expressed genes
    idx_zerosize <- apply(counts(sim_groups), MARGIN = 2, sum) == 0
    sim_groups <- sim_groups[, !idx_zerosize]

    if(scenario == 2){
      prop <- c(0.85, 0.05, 0.05,0.05)
    }else if(scenario == 16){
      prop <- c(0.45, 0.45, 0.05,0.05)
    }else if(scenario == 17){
      prop <- c(0.60, 0.30, 0.05,0.05)
    }else if(scenario == 18){
      prop <- c(0.35, 0.30, 0.30,0.05)
    }

  # 3.generate cell types
  print("generate cell types")

    ztrue <- label
    c_simu <- rep(NA, length(ztrue))
    if(scenario == 1){ # scenario 1: 4 cell types
      print("scenario == 1")
      c_simu <- ztrue # cell type assignment same as region assignment
    }else if(scenario != 1){
      print("scenario != 1")
      for(z in unique(label)){
        zi_idx <- ztrue == z # zi_idx is index for region z
        c_simu[zi_idx] <- sample(map_z2c(z), sum(zi_idx), prob = prop, replace = T)
      }
    }

    # 4.assign count data
    groups <- as.data.frame(colData(sim_groups)) %>% filter(Batch == "Batch" %&% 1)
    sim_cnt <- array(NA, c(J, N))
    for(c in 1:4){
      c_size <- sum(c_simu == c)  # true number of cells in cell type c
      c_cells <- groups$Cell[grepl(c, groups$Group)] # generated cells in group c
      cells_select <- sample(as.character(c_cells), c_size, replace = F)
      # sample same number of group c cells in real data from generated cells
      sim_cnt[, c_simu == c] <- as.matrix(counts(sim_groups)[, cells_select])
      # for positions of original cell type c cells, assign generated counts of group c
    }
    colnames(sim_cnt) <- "Cell" %&% 1:N
    rownames(sim_cnt) <- "Gene" %&% 1:J

  return(list(sim_cnt, c_simu, sim_seed))
}

map_z2c = function(z)
{
  case_when(
    z == 1 ~ c(1, 2, 3, 4),
    z == 2 ~ c(2, 3, 4, 1),
    z == 3 ~ c(3, 4, 1, 2),
    z == 4 ~ c(4, 1, 2, 3)
    )
}

```

















