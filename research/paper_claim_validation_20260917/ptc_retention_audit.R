#!/usr/bin/env Rscript
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
suppressPackageStartupMessages(library(jsonlite));suppressPackageStartupMessages(library(digest))
root<-'/fs/scratch/PCON0080/yimin/dgscrna'
base<-file.path(root,'results/hvg_ptc_20260916_v1/paper_claim_validation_20260917')
old<-file.path(root,'results/hvg_ptc_20260916_v1/r_reference_campaign_20260917')
gate<-fromJSON(file.path(base,'GBM_full_summary/GBM_FULL_DELIVERED.json'))
stopifnot(isTRUE(gate$PTC_compute_may_start),gate$receipt_sha256==digest(file=file.path(base,'GBM_full_summary/DELIVERY_RECEIPT.json'),algo='sha256'))
a<-commandArgs(trailingOnly=TRUE);cfg<-fromJSON(a[[1]]);route<-a[[2]]
stopifnot(cfg$kind=='retention')
source<-file.path(cfg$dest,route)
original<-file.path(old,'PTC_ablation',paste0('PTC_',cfg$group,'_GEOMETRY',cfg$budget,'_FIXED_CCAall_DL2000'),route)
fixed<-readLines(file.path(old,'PTC_geometry_fixed_inputs',cfg$group,'DL_features.txt'))
x<-readRDS(file.path(source,'DEG.rds'));y<-readRDS(file.path(original,'DEG.rds'))
x<-x[x$gene %in% fixed,,drop=FALSE]
key<-function(d)paste(as.character(d$cluster),d$gene,sep='\r')
stopifnot(!anyDuplicated(key(x)),!anyDuplicated(key(y)),setequal(key(x),key(y)))
x<-x[match(key(y),key(x)),,drop=FALSE]
columns<-intersect(c('p_val','avg_log2FC','pct.1','pct.2'),intersect(names(x),names(y)))
maxdiff<-list()
for(name in columns) {
  stopifnot(all(is.finite(x[[name]])),all(is.finite(y[[name]])))
  maximum<-if(nrow(y))max(abs(x[[name]]-y[[name]])) else 0
  stopifnot(maximum<=1e-12);maxdiff[[name]]<-maximum
}
write_json(list(status='passed',reference=original,shared_gene_DEG_rows=nrow(y),
  shared_gene_keys_exact=TRUE,unadjusted_statistics_max_absolute_difference=maxdiff,tolerance=1e-12,
  adjusted_p_values_not_compared='Gene-universe multiplicity changes; original scorer filters raw p and log2FC',
  original_DEG_sha256=digest(file=file.path(original,'DEG.rds'),algo='sha256'),
  new_DEG_sha256=digest(file=file.path(source,'DEG.rds'),algo='sha256'),job=Sys.getenv('SLURM_JOB_ID')),
  file.path(source,'retention_DEG_audit.json'),auto_unbox=TRUE,pretty=TRUE)
