#!/usr/bin/env Rscript
# Timeout recovery only: retain the installed FindAllMarkers implementation and
# restrict its outer identity loop; preserve every statistic, filter and order.
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
root <- '/fs/scratch/PCON0080/yimin/dgscrna'
out <- file.path(root,'results/hvg_ptc_20260916_v1/r_reference_campaign_20260917')
if(!nzchar(Sys.getenv('DGSCRNA_EXECUTION_SOURCE'))) {
  d<-file.path(out,'execution_sources',Sys.getenv('SLURM_JOB_ID'),Sys.getenv('SLURM_ARRAY_TASK_ID','single'))
  dir.create(d,recursive=TRUE,showWarnings=FALSE)
  p<-file.path(d,'checkpointed_deg.R')
  stopifnot(file.copy(file.path(root,'handoff/r_reference_campaign_20260917/checkpointed_deg.R'),p))
  Sys.setenv(DGSCRNA_EXECUTION_SOURCE=p);source(p,local=.GlobalEnv);quit(status=0)
}
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest))
options(future.globals.maxSize=32*1024^3)
plan(multicore,workers=4L)
args<-commandArgs(trailingOnly=TRUE)
mode<-args[[1]]
stopifnot(mode %in% c('validate','chunk','assemble'))
prep<-if(mode=='validate') file.path(out,'benchmark/brain_GBM/reference_CCA2000') else args[[2]]
route<-if(mode=='validate') 'PCA30_SNN' else args[[3]]
dest<-file.path(prep,route)
cl<-read.csv(file.path(dest,'clusters.csv'),colClasses='character')
identities<-factor(cl$cluster,levels=sort(unique(cl$cluster)))
ids<-as.character(sort(unique(identities)))
stopifnot(length(ids)>1L,!anyNA(identities),!anyDuplicated(cl$cell_id))
legacy_mean<-function(x) log2(Matrix::rowMeans(expm1(x))+1)
original<-Seurat::FindAllMarkers
original_hash<-digest(paste(deparse(original),collapse='\n'),algo='sha256',serialize=FALSE)

one_cluster<-function(object,id) {
  stopifnot(id %in% as.character(Idents(object)))
  fun<-original
  statements<-as.list(body(fun))
  target<-which(vapply(statements,function(x)identical(x,quote(genes.de <- list())),logical(1)))
  stopifnot(length(target)==1L)
  # Insert after identity discovery, before the unchanged loops and filters.
  statements<-append(statements,list(quote({
    idents.all<-idents.all[as.character(idents.all) %in% selected_cluster_ids]
    stopifnot(length(idents.all)==1L)
  })),after=target-1L)
  body(fun)<-as.call(statements)
  environment(fun)<-list2env(list(selected_cluster_ids=id),parent=environment(original))
  fun(object,assay=DefaultAssay(object),slot='data',test.use='wilcox_limma',
    logfc.threshold=.25,min.pct=.1,min.diff.pct=-Inf,only.pos=FALSE,
    max.cells.per.ident=Inf,random.seed=1,min.cells.feature=3,min.cells.group=3,
    pseudocount.use=1,mean.fxn=legacy_mean,fc.name='avg_log2FC',base=2,
    return.thresh=.01,densify=FALSE,verbose=FALSE)
}
combine<-function(parts) {
  d<-data.frame()
  for(part in parts)if(nrow(part)>0L)d<-rbind(d,part)
  if(nrow(d)==0L)return(data.frame(cluster=character(),gene=character(),avg_log2FC=numeric()))
  rownames(d)<-make.unique(as.character(d$gene))
  stopifnot(all(is.finite(d$avg_log2FC)))
  d
}
atomic_rds<-function(x,p) {
  temporary<-paste0(p,'.part.',Sys.getenv('SLURM_JOB_ID'),'.',Sys.getenv('SLURM_ARRAY_TASK_ID','single'))
  saveRDS(x,temporary);stopifnot(file.rename(temporary,p))
}
atomic_json<-function(x,p) {
  temporary<-paste0(p,'.part.',Sys.getenv('SLURM_JOB_ID'),'.',Sys.getenv('SLURM_ARRAY_TASK_ID','single'))
  write_json(x,temporary,auto_unbox=TRUE,pretty=TRUE);stopifnot(file.rename(temporary,p))
}
source_hash<-digest(file=Sys.getenv('DGSCRNA_EXECUTION_SOURCE'),algo='sha256')
audit_path<-file.path(out,'verification/checkpointed_DEG_equivalence.json')

if(mode %in% c('validate','chunk')) {
  obj<-readRDS(file.path(prep,'expression_PCA30.rds'))
  stopifnot(identical(cl$cell_id,colnames(obj)))
  Idents(obj)<-identities
}
if(mode=='validate') {
  parts<-lapply(ids,function(id){
    cat('VALIDATE_CLUSTER',id,'\n');flush.console();one_cluster(obj,id)
  })
  observed<-combine(parts)
  reference<-readRDS(file.path(dest,'DEG.rds'))
  stopifnot(identical(observed,reference))
  atomic_json(list(status='passed',job=Sys.getenv('SLURM_JOB_ID'),
    rows=nrow(observed),clusters=length(ids),exact_all_fields_attributes_and_row_order=TRUE,
    installed_FindAllMarkers_sha256=original_hash,source_sha256=source_hash,
    reference_DEG_sha256=digest(file=file.path(dest,'DEG.rds'),algo='sha256'),
    scope='Every saved GBM PCA-SNN DEG field, factor attribute and row order equals the installed original full-loop output; only the outer cluster iteration is partitioned.'),audit_path)
  cat('EXACT_CHECKPOINTED_DEG_PARITY',nrow(observed),'\n');quit(status=0)
}

validation<-fromJSON(audit_path)
stopifnot(validation$status=='passed',validation$source_sha256==source_hash,
  validation$installed_FindAllMarkers_sha256==original_hash)
checkpoints<-file.path(dest,'DEG_checkpoints');dir.create(checkpoints,showWarnings=FALSE)
input_hashes<-list(expression=digest(file=file.path(prep,'expression_PCA30.rds'),algo='sha256'),
  clusters=digest(file=file.path(dest,'clusters.csv'),algo='sha256'))
cluster_path<-function(id)file.path(checkpoints,paste0('cluster_',digest(id,algo='sha256',serialize=FALSE),'.rds'))
read_part<-function(id) {
  p<-cluster_path(id);meta<-fromJSON(paste0(p,'.json'))
  stopifnot(identical(meta$cluster,id),identical(meta$input_hashes,input_hashes),
    meta$source_sha256==source_hash,meta$installed_FindAllMarkers_sha256==original_hash,
    meta$sha256==digest(file=p,algo='sha256'))
  d<-readRDS(p)
  stopifnot(nrow(d)==meta$rows,nrow(d)==0L || all(as.character(d$cluster)==id))
  d
}
if(mode=='chunk') {
  size<-as.integer(args[[4]])
  index<-as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
  stopifnot(size>=1L,index>=0L)
  selected<-ids[seq_along(ids)>index*size & seq_along(ids)<=(index+1L)*size]
  stopifnot(length(selected)>0L)
  for(id in selected) {
    p<-cluster_path(id)
    if(file.exists(p) && file.exists(paste0(p,'.json'))) {
      read_part(id);cat('REUSED_CLUSTER',id,'\n');flush.console();next
    }
    cat(format(Sys.time()),'START_CLUSTER',id,'\n');flush.console()
    d<-one_cluster(obj,id)
    atomic_rds(d,p)
    atomic_json(list(cluster=id,rows=nrow(d),input_hashes=input_hashes,
      source_sha256=source_hash,installed_FindAllMarkers_sha256=original_hash,
      sha256=digest(file=p,algo='sha256'),job=Sys.getenv('SLURM_JOB_ID'),
      task=Sys.getenv('SLURM_ARRAY_TASK_ID')),paste0(p,'.json'))
    cat(format(Sys.time()),'COMPLETE_CLUSTER',id,nrow(d),'\n');flush.console()
  }
} else {
  d<-combine(lapply(ids,read_part))
  target<-file.path(dest,'DEG.rds')
  if(file.exists(target))stopifnot(identical(readRDS(target),d)) else atomic_rds(d,target)
  atomic_json(list(status='assembled',job=Sys.getenv('SLURM_JOB_ID'),clusters=length(ids),
    rows=nrow(d),input_hashes=input_hashes,source_sha256=source_hash,
    installed_FindAllMarkers_sha256=original_hash,sha256=digest(file=target,algo='sha256'),
    cluster_order=ids,validation_sha256=digest(file=audit_path,algo='sha256')),
    file.path(checkpoints,'ASSEMBLED.json'))
  cat('ASSEMBLED_DEG',length(ids),nrow(d),'\n')
}
