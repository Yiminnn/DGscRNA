stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
root<-'/fs/scratch/PCON0080/yimin/dgscrna'
out<-file.path(root,'results/hvg_ptc_20260916_v1/r_reference_campaign_20260917')
if(!nzchar(Sys.getenv('DGSCRNA_EXECUTION_SOURCE'))) {
  d<-file.path(out,'execution_sources',Sys.getenv('SLURM_JOB_ID'));dir.create(d,recursive=TRUE,showWarnings=FALSE)
  p<-file.path(d,'audit_score_recovery_and_workers.R')
  stopifnot(file.copy(file.path(root,'handoff/r_reference_campaign_20260917/audit_score_recovery_and_workers.R'),p))
  Sys.setenv(DGSCRNA_EXECUTION_SOURCE=p);source(p,local=.GlobalEnv);quit(status=0)
}
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest))
driver<-file.path(root,'handoff/r_reference_campaign_20260917/reference_score.R')
for(s in parse(driver))if(is.call(s)&&identical(s[[1]],as.name('<-'))&&
  as.character(s[[2]]) %in% c('publish','csv_atomic','rds_atomic','read_cache'))eval(s,envir=.GlobalEnv)
prep<-file.path(out,'verification/score_cache_recovery_test',Sys.getenv('SLURM_JOB_ID'));dir.create(prep,recursive=TRUE,showWarnings=FALSE)
tab<-data.frame(cell_id=c('a','b','c'),cluster=c('0','1','1'))
csv<-file.path(prep,'clusters.csv');csv_atomic(tab,csv,row.names=FALSE)
valid<-function(x)identical(x,tab)
stopifnot(identical(read_cache(csv,function(p)read.csv(p,colClasses='character'),valid),tab))
writeLines(c('"cell_id","cluster"','"a","0"'),csv)
stopifnot(is.null(read_cache(csv,function(p)read.csv(p,colClasses='character'),valid)),!file.exists(csv))
csv_atomic(tab,csv,row.names=FALSE)
rds<-file.path(prep,'DEG.rds');rds_atomic(tab,rds)
original_bytes<-readBin(rds,'raw',n=file.info(rds)$size)
writeBin(original_bytes[seq_len(17)],rds)
stopifnot(is.null(read_cache(rds,readRDS,valid)),!file.exists(rds))
rds_atomic(tab,rds);stopifnot(identical(read_cache(rds,readRDS,valid),tab))
stopifnot(length(list.files(prep,pattern='invalid_'))==2L,
  length(list.files(prep,pattern='\\.part\\.'))==0L)

case<-file.path(out,'benchmark/brain_GBM/reference_CCA2000')
obj<-readRDS(file.path(case,'expression_PCA30.rds'))
cl<-read.csv(file.path(case,'PCA30_SNN/clusters.csv'),colClasses='character')
stopifnot(identical(cl$cell_id,colnames(obj)))
Idents(obj)<-factor(cl$cluster,levels=sort(unique(cl$cluster)))
legacy_mean<-function(x)log2(Matrix::rowMeans(expm1(x))+1)
normalise<-function(d){d$cluster<-as.character(d$cluster);d<-d[order(d$cluster,d$gene),,drop=FALSE];rownames(d)<-NULL;d}
saved<-normalise(readRDS(file.path(case,'PCA30_SNN/DEG.rds')))
checks<-list();warnings_seen<-list();options(future.globals.maxSize=16*1024^3)
for(workers in c(1L,2L,4L)) {
  if(workers==1L)plan(sequential) else plan(multicore,workers=workers)
  warnings_this<-character()
  d<-withCallingHandlers(FindAllMarkers(obj,assay=DefaultAssay(obj),slot='data',test.use='wilcox_limma',
    logfc.threshold=.25,min.pct=.1,min.diff.pct=-Inf,only.pos=FALSE,
    max.cells.per.ident=Inf,random.seed=1,min.cells.feature=3,min.cells.group=3,
    pseudocount.use=1,mean.fxn=legacy_mean,fc.name='avg_log2FC',base=2,
    return.thresh=.01,densify=FALSE,verbose=FALSE),warning=function(w){
      warnings_this<<-c(warnings_this,conditionMessage(w));invokeRestart('muffleWarning')})
  d<-normalise(d)
  comparison<-all.equal(saved,d,tolerance=0,check.attributes=TRUE)
  stopifnot(isTRUE(comparison))
  checks[[as.character(workers)]]<-list(workers=workers,DEG_rows=nrow(d),exact_all_fields_to_saved=TRUE)
  warnings_seen[[as.character(workers)]]<-unique(warnings_this)
  cat('EXACT_DEG_PARITY',workers,nrow(d),'\n');flush.console()
}
plan(sequential)
write_json(list(status='passed',job=Sys.getenv('SLURM_JOB_ID'),
  cache_recovery=list(truncated_CSV_preserved_and_recomputed=TRUE,truncated_RDS_preserved_and_recomputed=TRUE,atomic_publication_verified=TRUE),
  worker_equivalence=checks,warnings=warnings_seen,
  scope='Full saved GBM PCA-SNN expression and partition; deterministic Wilcoxon/mean/filter parameters identical to campaign. One, two and four workers reproduce every saved DEG field exactly.',
  source_sha256=digest(file=Sys.getenv('DGSCRNA_EXECUTION_SOURCE'),algo='sha256'),
  scorer_sha256=digest(file=driver,algo='sha256')),
  file.path(out,'verification/score_recovery_and_DEG_workers.json'),auto_unbox=TRUE,pretty=TRUE)
