#!/usr/bin/env Rscript
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest))
args<-commandArgs(trailingOnly=TRUE)
cfg<-fromJSON(args[[1]])
prep<-cfg$prep; dest<-cfg$geometry
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
hash<-function(p)digest(file=p,algo='sha256')
pm<-fromJSON(file.path(prep,'prepare_manifest.json'))
stopifnot(readLines(file.path(prep,'PREPARED'))==hash(file.path(prep,'prepare_manifest.json')))
stopifnot(hash(file.path(prep,'expression_PCA30.rds'))==pm$expression_sha256)
obj<-readRDS(file.path(prep,'expression_PCA30.rds'))
cells<-read.csv(file.path(prep,'cells.csv'),stringsAsFactors=FALSE)
stopifnot(identical(colnames(obj),cells$cell_id),DefaultAssay(obj)=='RNA')
x<-obj[['RNA']]@scale.data
stopifnot(identical(colnames(x),cells$cell_id),all(is.finite(x)),
          setequal(rownames(x),readLines(file.path(prep,'selected_features.txt'))))
target<-file.path(dest,'scaled_HVG.float64.bin')
con<-file(paste0(target,'.part'),'wb')
# R gene-by-cell column-major order is Python cell-by-gene row-major order.
writeBin(as.numeric(x),con,size=8,endian='little');close(con)
stopifnot(file.rename(paste0(target,'.part'),target))
write.csv(cells,file.path(dest,'cells.csv'),row.names=FALSE)
writeLines(rownames(x),file.path(dest,'features.txt'))
script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[[1]])
m<-list(status='completed',n_cells=ncol(x),n_features=nrow(x),binary=target,
  binary_sha256=hash(target),cells_sha256=hash(file.path(dest,'cells.csv')),
  features_sha256=hash(file.path(dest,'features.txt')),
  prepare_manifest_sha256=hash(file.path(prep,'prepare_manifest.json')),
  original_expression_sha256=pm$expression_sha256,
  input='Exact native R ScaleData selected-HVG matrix; centered/scaled; scale.max=10',
  dtype='little-endian float64; cell-major',source_sha256=hash(script),
  job=Sys.getenv('SLURM_JOB_ID'),reference_labels_used_for_fit=FALSE)
write_json(m,file.path(dest,'manifest.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(hash(file.path(dest,'manifest.json')),file.path(dest,'COMPLETE'))
