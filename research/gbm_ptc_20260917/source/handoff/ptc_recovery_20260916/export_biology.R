#!/usr/bin/env Rscript
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest))
base <- '/fs/scratch/PCON0080/yimin/dgscrna/results/hvg_ptc_20260916_v1/ptc_experiments'
group <- c('MTN','TUT')[[as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))+1L]]
stopifnot(file.exists(file.path(base,'protocol/LABEL_RULES_FROZEN')))
modules <- fromJSON(file.path(base,'protocol/RNA_modules.json'))
obj <- readRDS(file.path(base,'prepared',group,'RNA_normalized.rds'))
dest <- file.path(base,'biology');dir.create(dest,showWarnings=FALSE)
out <- data.frame(cell_id=colnames(obj),stringsAsFactors=FALSE)
present <- list()
for(name in names(modules)) {
  genes <- intersect(modules[[name]],rownames(obj));stopifnot(length(genes)>0)
  present[[name]] <- genes
  out[[paste0(name,'_mean_log1p')]] <- Matrix::colMeans(obj[['RNA']]@data[genes,,drop=FALSE])
  out[[paste0(name,'_n_detected')]] <- Matrix::colSums(obj[['RNA']]@counts[genes,,drop=FALSE]>0)
}
out$T_RNA_support_ge2 <- out$T_core_n_detected>=2L
genes <- unique(unlist(modules,use.names=FALSE));genes <- intersect(genes,rownames(obj))
expr <- t(as.matrix(obj[['RNA']]@data[genes,,drop=FALSE]))
con <- gzfile(file.path(dest,paste0(group,'.gene_expression.csv.gz')),'wt')
write.csv(data.frame(cell_id=rownames(expr),expr,check.names=FALSE),con,row.names=FALSE);close(con)
con <- gzfile(file.path(dest,paste0(group,'.module_scores.csv.gz')),'wt');write.csv(out,con,row.names=FALSE);close(con)
write_json(list(group=group,n_cells=ncol(obj),RNA_genes=nrow(obj),present_module_genes=present,
   score='Mean unscaled RNA LogNormalize10k log1p expression and number of genes with nonzero raw counts',
   independent_truth=FALSE,source_sha256=digest(file=Sys.getenv('PTC_BIO_SOURCE'),algo='sha256'),
   modules_sha256=digest(file=file.path(base,'protocol/RNA_modules.json'),algo='sha256'),job=Sys.getenv('SLURM_JOB_ID')),
   file.path(dest,paste0(group,'.manifest.json')),pretty=TRUE,auto_unbox=TRUE)
writeLines(digest(file=file.path(dest,paste0(group,'.manifest.json')),algo='sha256'),file.path(dest,paste0(group,'.COMPLETE')))
