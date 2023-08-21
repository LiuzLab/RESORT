# setwd("/mnt/atlas_local/linhua/data/ReSort_Revise/simulation/scripts/run_deconv_models/")
# install the MuSiC package
# install.packages('devtools')
library(devtools)
devtools::install_github('xuranw/MuSiC')
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biobase")
library(MuSiC)
library(ggplot2)
library(ggpubr)
library(Biobase)

runMusic <- function(mixture_count_fn, ref_count_fn, ref_meta_fn, col){
  mixture = read.csv(mixture_count_fn, row.names=1)
  mixture = t(mixture)
  mixture.eset <- ExpressionSet(assayData = mixture)
  
  ref_count <-  t(as.matrix(read.csv(ref_count_fn, row.names=1, 
                                      check.names = F, stringsAsFactors = F)))
  ref_meta <- read.csv(ref_meta_fn, row.names=1, 
                        check.names = F, stringsAsFactors = F)
  
  
  ref_meta$spotID <- row.names(ref_meta)
  ref_meta <- ref_meta[, c(col, 'spotID')]
  ref_meta.meta <- data.frame(labelDescription=c(col, "spotID"),
                               row.names=c(col, "spotID"))
  
  phenoData <- new("AnnotatedDataFrame",data=ref_meta, varMetadata=ref_meta.meta)
  
  ref.eset <- ExpressionSet(assayData = ref_count, phenoData = phenoData)
  Est.prop.ref = music_prop(bulk.eset = mixture.eset, sc.eset = ref.eset,
                             clusters = col, samples = 'spotID',
                             select.ct = NULL, 
                             verbose = T)
  props <- as.data.frame(Est.prop.ref$Est.prop.weighted)
  return (props)
}

# ### INPUT FILES, TO BE CHANGED ####
# projDir <- '/houston_20t/alexw/ST/iSort/simulation_finer_resolution/data/cancer_inf_0.1/'
# outDir <- file.path(projDir, 'music_minor_results/')
# workDir <- getwd()
# ref_count_fn <- file.path(workDir, 'single_cell_ref/single_cell_ref_count_corrected_augmented_external.csv')
# ref_meta_fn <- file.path(workDir, 'single_cell_ref/single_cell_ref_meta_corrected_augmented_external.csv')
# mixture_count_fn <- file.path(projDir, 'simulated_mixture_raw_counts.csv')
# mixture.props <- runMusic(mixture_count_fn, ref_count_fn, ref_meta_fn, 'bio_celltype')
# write.csv(mixture.props, file.path(projDir, 'music_minor_results/estimated_proportions_external_corrected_augmented.csv'))
# # system.time(runMusic_project(projDir, outDir, col = 'bio_celltype'))
# # system.time(runMusic_project(projDir, col = 'bio_celltype_major'))
