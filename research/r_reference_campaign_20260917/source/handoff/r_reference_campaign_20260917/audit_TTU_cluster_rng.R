#!/usr/bin/env Rscript
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
root<-'/fs/scratch/PCON0080/yimin/dgscrna'
out<-file.path(root,'results/hvg_ptc_20260916_v1/r_reference_campaign_20260917')
if(!nzchar(Sys.getenv('DGSCRNA_EXECUTION_SOURCE'))) {
  d<-file.path(out,'execution_sources',Sys.getenv('SLURM_JOB_ID'));dir.create(d,recursive=TRUE,showWarnings=FALSE)
  p<-file.path(d,'audit_TTU_cluster_rng.R')
  stopifnot(file.copy(file.path(root,'handoff/r_reference_campaign_20260917/audit_TTU_cluster_rng.R'),p))
  Sys.setenv(DGSCRNA_EXECUTION_SOURCE=p);source(p,local=.GlobalEnv);quit(status=0)
}
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest))
options(future.globals.maxSize=32*1024^3)
prep<-file.path(out,'PTC_ablation/PTC_TTU_CCAall')
base<-readRDS(file.path(prep,'expression_PCA30.rds'))
ref<-read.csv(file.path(prep,'UMAP2_SNN/clusters.csv'),colClasses='character')
z<-as.matrix(read.csv(file.path(prep,'UMAP2.csv'),row.names=1,check.names=FALSE))
stopifnot(identical(ref$cell_id,colnames(base)),identical(rownames(z),colnames(base)))
colnames(z)<-paste0('UMAP_',1:2)
base[['umap']]<-CreateDimReducObject(embeddings=z,key='UMAP_',assay=DefaultAssay(base))
canonical<-function(x)match(x,unique(x))
checks<-list();warnings_seen<-list()
capture_phase<-function(phase,expr)withCallingHandlers(expr,warning=function(w) {
  warnings_seen[[length(warnings_seen)+1L]]<<-list(phase=phase,workers=workers,
    global_seed=seed_name,message=conditionMessage(w))
})
for(workers in c(1L,4L))for(seed in c(NA_integer_,42L,2026L)) {
  if(workers==1L)plan(sequential) else plan(multicore,workers=workers)
  seed_name<-if(is.na(seed))'unset' else as.character(seed)
  if(is.na(seed)) {
    if(exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE))rm('.Random.seed',envir=.GlobalEnv)
  } else set.seed(seed)
  obj<-capture_phase('FindNeighbors',FindNeighbors(base,reduction='umap',dims=1:2,k.param=20L,
    compute.SNN=TRUE,prune.SNN=1/15,nn.method='annoy',n.trees=50,annoy.metric='euclidean',
    graph.name=c('ref_nn','ref_snn'),verbose=FALSE))
  obj<-capture_phase('FindClusters',FindClusters(obj,graph.name='ref_snn',resolution=.5,
    algorithm=1,modularity.fxn=1,n.start=10,n.iter=10,random.seed=0,verbose=FALSE))
  observed<-as.character(Idents(obj))
  exact<-identical(observed,ref$cluster)
  equivalent<-identical(canonical(observed),canonical(ref$cluster))
  key<-paste(workers,seed_name,sep='_')
  checks[[key]]<-list(workers=workers,global_seed=seed_name,Louvain_seed=0,
    raw_labels_identical=exact,partition_identical_up_to_label_ids=equivalent,
    raw_label_mismatches=sum(observed!=ref$cluster),clusters=length(unique(observed)))
  if(!equivalent)write.csv(data.frame(cell_id=ref$cell_id,reference=ref$cluster,replayed=observed),
    file.path(out,'verification',paste0('TTU_cluster_rng_mismatch_',key,'.csv')),row.names=FALSE)
  rm(obj);invisible(gc())
  cat('TTU_CLUSTER_REPLAY',key,'exact_labels',exact,'exact_partition',equivalent,'\n');flush.console()
}
plan(sequential)
passed<-all(vapply(checks,function(x)x$partition_identical_up_to_label_ids,logical(1)))
phases<-unique(vapply(Filter(function(x)grepl('UNRELIABLE VALUE',x$message),warnings_seen),function(x)x$phase,character(1)))
source_text<-unlist(lapply(c('FindClusters.default','RunModularityClustering','GroupSingletons'),function(name)
  c(name,deparse(get(name,asNamespace('Seurat'))))))
writeLines(source_text,file.path(out,'verification/TTU_clustering_installed_source.txt'))
result<-list(status=if(passed)'passed' else 'failed',job=Sys.getenv('SLURM_JOB_ID'),
  observed_warning_job='7351708',cells=ncol(base),checks=checks,warnings=warnings_seen,
  replicated_future_RNG_warning_phases=phases,
  input_clusters_sha256=digest(file=file.path(prep,'UMAP2_SNN/clusters.csv'),algo='sha256'),
  UMAP_sha256=digest(file=file.path(prep,'UMAP2.csv'),algo='sha256'),
  source_sha256=digest(file=Sys.getenv('DGSCRNA_EXECUTION_SOURCE'),algo='sha256'),
  scope='Original saved TTU all-gene UMAP coordinates, identical SNN and Louvain parameters. One/four workers; global RNG unset/42/2026; original Louvain random.seed=0 unchanged. Outputs are compared to the saved actual partition; nothing overwrites it.')
temporary<-file.path(out,'verification/TTU_cluster_rng_audit.json.part')
write_json(result,temporary,auto_unbox=TRUE,pretty=TRUE)
stopifnot(file.rename(temporary,file.path(out,'verification/TTU_cluster_rng_audit.json')),passed)
