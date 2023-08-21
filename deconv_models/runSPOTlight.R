devtools::install_github("https://github.com/MarcElosua/SPOTlight")
#install.packages("dplyr")
#install.packages('Seurat')
#workDir <- "/houston_20t/alexw/ST/iSort/simulation_finer_resolution/"
setwd(workDir)

library(Seurat)
library(SPOTlight)
library(dplyr)

run_spotlight <- function(mixture_fn, ref_count_fn, ref_meta_fn, col) {
  mixtures <- t(as.matrix(read.csv(mixture_fn,row.names=1)))
  ## PROCESS REFERENCE DATA
  ref <-  t(as.matrix(read.csv(ref_count_fn, row.names=1, 
                               check.names = F, stringsAsFactors = F)))
  ref_meta <- read.csv(ref_meta_fn, row.names=1, 
                       check.names = F, stringsAsFactors = F)
  cts <- ref_meta[, col]
  colnames(ref) <- paste("spot", 1:dim(ref)[2], sep="")
  pheno_sigs <- data.frame(row.names=colnames(ref), spotID = colnames(ref), class=cts)
  sc_obj <- CreateSeuratObject(count=ref, project='sc_ref', meta.data = pheno_sigs)
  sc_obj <- NormalizeData(sc_obj, normalization.method = "LogNormalize",
                          scale.factor = 1e6, margin = 1)
  Seurat::Idents(object = sc_obj) <- sc_obj@meta.data$class
  cluster_markers_all <- Seurat::FindAllMarkers(object = sc_obj, assay = "RNA",
                                                slot = "data", verbose = TRUE, only.pos = TRUE)
  cluster_markers_all <- cluster_markers_all[cluster_markers_all$avg_log2FC >= 0.26 &
                                               cluster_markers_all$p_val_adj <= 0.05, ]
  ## PROCESS SPATIAL DATA
  spatial_counts <- CreateSeuratObject(count=mixtures, project='st_mixture')
  spatial_counts <- NormalizeData(spatial_counts, normalization.method = 'LogNormalize', 
                                  scale.factor = 1e6, margin=1)
  
  ## RUN SPOTLIGHT
  spotlight_ls <- spotlight_deconvolution(
    se_sc = sc_obj,
    counts_spatial = spatial_counts@assays$RNA@counts,
    clust_vr = "class", # Variable in sc_seu containing the cell-type annotation
    cluster_markers = cluster_markers_all, # Dataframe with the marker genes
    ntop = NULL, # How many of the marker genes to use (by default all)
    transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
    method = "nsNMF", # Factorization method
    min_cont = 0 # Remove those cells contributing to a spot below a certain threshold 
  )
  
  ## Post-process results
  spotlight.props <- data.frame(spotlight_ls[[2]], row.names=colnames(mixtures))
  xs <- c()
  ys <- c()
  for (spot in row.names(spotlight.props)){
    x <- as.integer(strsplit(spot, split = "x")[[1]][1])
    y <- as.integer(strsplit(spot, split = "x")[[1]][2])
    xs <- c(xs, x)
    ys <- c(ys, y)
  }
  
  spotlight.props$coordX <- as.numeric(xs)
  spotlight.props$coordY <- as.numeric(ys)
  spotlight.props <- spotlight.props[, c(unique(cts), c('coordX', 'coordY'))]
  return(spotlight.props)
}