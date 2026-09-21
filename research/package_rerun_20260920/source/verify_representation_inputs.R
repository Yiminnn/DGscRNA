stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
extra<-Sys.getenv('DGSCRNA_REFERENCE_R_LIB','')
.libPaths(c(if(nzchar(extra))extra else character(),.Library),include.site=FALSE)
suppressPackageStartupMessages(library(Seurat))
a<-commandArgs(trailingOnly=TRUE)
stopifnot(length(a)==5L)
actual<-readRDS(file.path(a[[1]],'expression_PCA30.rds'))
source<-readRDS(file.path(a[[2]],'expression_PCA30.rds'))
expected<-readRDS(file.path(a[[3]],'expression_PCA30.rds'))
stopifnot(identical(actual[['RNA']],source[['RNA']]),
          identical(actual[['RNA']],expected[['RNA']]),
          identical(colnames(actual),colnames(source)),
          identical(colnames(actual),colnames(expected)))
Z<-Embeddings(actual,'control')
if(a[[5]]=='default') {
  reduction<-if(a[[4]]=='PCA30')'pca' else 'umap'
  stopifnot(identical(unname(Z),unname(Embeddings(source,reduction))),
            identical(unname(Z),unname(Embeddings(expected,reduction))))
} else stopifnot(identical(Z,Embeddings(expected,'control')))
genes<-if(a[[4]]=='RNA_noDR')rownames(actual[['RNA']]@scale.data) else rownames(Loadings(actual,'pca'))
stopifnot(identical(readLines(file.path(a[[1]],'control_feature_genes.txt')),genes))
for(name in c('selected_features.txt','geometry_features.txt','scoring_features.txt','DL_features.txt'))
  stopifnot(identical(readLines(file.path(a[[1]],name)),readLines(file.path(a[[2]],name))),
            identical(readLines(file.path(a[[1]],name)),readLines(file.path(a[[3]],name))))
cat('RNA_ASSAY_AND_GENE_ORDER_EXACT_CONTROL_COORDINATES_EXACT\n')
