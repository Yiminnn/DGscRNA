stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(jsonlite))
root<-'/fs/scratch/PCON0080/yimin/dgscrna'
base<-file.path(root,'results/hvg_ptc_20260916_v1/r_reference_campaign_20260917')
prep<-file.path(base,'benchmark/brain_GBM/reference_CCA2000')
obj<-readRDS(file.path(prep,'expression_PCA30.rds'))
route<-file.path(prep,'PCA30_SNN')
cl<-read.csv(file.path(route,'clusters.csv'),stringsAsFactors=FALSE)
stopifnot(identical(cl$cell_id,colnames(obj)))
obj$reference_clusters<-as.character(cl$cluster)
libs<-fromJSON(file.path(base,'markers/brain_GBM.json'),simplifyVector=FALSE)
libs<-lapply(libs,function(l)lapply(l,unlist,use.names=FALSE))
deg<-readRDS(file.path(route,'DEG.rds'))
initial<-read.csv(gzfile(file.path(route,'initial_calls.csv.gz')),stringsAsFactors=FALSE,check.names=FALSE)
source<-file.path(root,'results/hvg_ptc_20260916_v1/ptc_recovery/archive/tcr/ptc_val/scripts/DGscRNA-Share/R/source.R')
for(statement in parse(source)) if(is.call(statement) && identical(statement[[1]],as.name('<-')) &&
  identical(statement[[2]],as.name('density_score')))eval(statement,envir=.GlobalEnv)
checks<-list()
for(i in seq_along(libs)) {
  lib<-names(libs)[[i]]
  obj<-density_score(obj,markers=libs[[i]],DEG_markers_set=list(reference_clusters=deg),
      annotation_name=lib,clusterings='reference_clusters',cutoffs=c('none','mean','0.5'))
  for(cut in c('none','mean','0.5')) {
    aid<-sprintf('L%02d_%s',i-1L,if(cut=='0.5')'p050' else cut)
    got<-obj@meta.data[[paste(lib,'reference_clusters',cut,sep='_')]]
    stopifnot(identical(as.character(got),as.character(initial[[aid]])))
    checks[[aid]]<-list(library=lib,cutoff=cut,n_cells=ncol(obj),exact=TRUE)
  }
}
write_json(list(status='passed',checks=checks,job=Sys.getenv('SLURM_JOB_ID'),
 criterion='Every GBM PCA-SNN marker/cutoff initial call equals unmodified archived density_score'),
 file.path(base,'density_original_function_parity.json'),auto_unbox=TRUE,pretty=TRUE)
cat('All 24 density conditions exactly match original function\n')
