stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest))
root<-'/fs/scratch/PCON0080/yimin/dgscrna'
base<-file.path(root,'results/hvg_ptc_20260916_v1/r_reference_campaign_20260917')
unit<-readLines(file.path(base,'pre_guard_audit_units.txt'))[[as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))+1L]]
prep<-file.path(base,'benchmark',unit,'reference_CCA2000')
original<-file.path(root,'results/hvg_ptc_20260916_v1/ptc_recovery/archive/tcr/ptc_val/scripts/DGscRNA-Share/R/source.R')
for(statement in parse(original))if(is.call(statement)&&identical(statement[[1]],as.name('<-'))&&identical(statement[[2]],as.name('density_score')))eval(statement,envir=.GlobalEnv)
libs<-fromJSON(file.path(base,'markers',paste0(unit,'.json')),simplifyVector=FALSE)
libs<-lapply(libs,function(l)lapply(l,unlist,use.names=FALSE))
checks<-list()
for(route in c('PCA30_SNN','PCA30_HDBSCAN_R','UMAP2_SNN','UMAP2_HDBSCAN_R')) {
  source<-file.path(prep,route)
  cl<-read.csv(file.path(source,'clusters.csv'),stringsAsFactors=FALSE,colClasses='character')
  ids<-sort(unique(cl$cluster));index<-setNames(seq_along(ids)-1L,ids)
  # The historical function assumes contiguous zero-based IDs. Relabel only
  # those arbitrary cluster IDs for this audit; no cells, DEGs or type labels change.
  deg<-readRDS(file.path(source,'DEG.rds'));deg$cluster<-unname(index[as.character(deg$cluster)])
  initial<-read.csv(gzfile(file.path(source,'initial_calls.csv.gz')),stringsAsFactors=FALSE,check.names=FALSE,colClasses='character')
  stopifnot(identical(initial$cell_id,cl$cell_id))
  dummy<-matrix(1,nrow=2,ncol=length(ids),dimnames=list(c('AUDIT_A','AUDIT_B'),paste0('cluster_',ids)))
  obj<-CreateSeuratObject(dummy,min.cells=0,min.features=0)
  obj$reference_clusters<-as.character(seq_along(ids)-1L)
  for(i in seq_along(libs)) {
    lib<-names(libs)[[i]]
    obj<-density_score(obj,markers=libs[[i]],DEG_markers_set=list(reference_clusters=deg),annotation_name=lib,
                        clusterings='reference_clusters',cutoffs=c('none','mean','0.5'))
    for(cut in c('none','mean','0.5')) {
      aid<-sprintf('L%02d_%s',i-1L,if(cut=='0.5')'p050' else cut)
      calls<-obj@meta.data[[paste(lib,'reference_clusters',cut,sep='_')]]
      got<-as.character(calls[unname(index[cl$cluster])+1L])
      stopifnot(identical(got,as.character(initial[[aid]])))
      checks[[paste(route,aid,sep='/')]]<-list(exact=TRUE,n_cells=nrow(cl),n_clusters=length(ids),library=lib,cutoff=cut)
    }
  }
  cat('PARITY',unit,route,length(libs)*3L,'\n');flush.console()
}
dest<-file.path(base,'verification/pre_guard_density');dir.create(dest,recursive=TRUE,showWarnings=FALSE)
write_json(list(status='passed',unit=unit,job=Sys.getenv('SLURM_JOB_ID'),checks=checks,
   original_function_source_sha256=digest(file=original,algo='sha256'),
   meaning='Independent unmodified original density_score reproduces every initial call on the saved DEG/partition inputs.',
   provenance_limit='Validates outputs; does not reconstruct the exact transient script bytes of pre-guard jobs.'),
   file.path(dest,paste0(unit,'.json')),auto_unbox=TRUE,pretty=TRUE)
