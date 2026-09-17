#!/usr/bin/env Rscript
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
root<-'/fs/scratch/PCON0080/yimin/dgscrna'
out<-file.path(root,'results/hvg_ptc_20260916_v1/r_reference_campaign_20260917')
if(!nzchar(Sys.getenv('DGSCRNA_EXECUTION_SOURCE'))) {
  d<-file.path(out,'execution_sources',Sys.getenv('SLURM_JOB_ID'));dir.create(d,recursive=TRUE,showWarnings=FALSE)
  p<-file.path(d,'audit_TTU_wilcox_kernel.R')
  stopifnot(file.copy(file.path(root,'handoff/r_reference_campaign_20260917/audit_TTU_wilcox_kernel.R'),p))
  Sys.setenv(DGSCRNA_EXECUTION_SOURCE=p);source(p,local=.GlobalEnv);quit(status=0)
}
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(future))
suppressPackageStartupMessages(library(future.apply))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest))
prep<-file.path(out,'PTC_ablation/PTC_TTU_CCAall')
dest<-file.path(prep,'UMAP2_SNN')
deg<-readRDS(file.path(dest,'DEG.rds'))
cl<-read.csv(file.path(dest,'clusters.csv'),colClasses='character')
obj<-readRDS(file.path(prep,'expression_PCA30.rds'))
stopifnot(identical(cl$cell_id,colnames(obj)))
ncells<-ncol(obj)
stopifnot(!is.na(suppressWarnings(ncells*ncells))) # The actual Seurat limma branch.
X<-GetAssayData(obj,assay=DefaultAssay(obj),layer='data')
ids<-sort(unique(as.character(deg$cluster)))
cases<-do.call(rbind,lapply(ids,function(id) {
  d<-deg[as.character(deg$cluster)==id,,drop=FALSE]
  chosen<-unique(c(1L,ceiling(nrow(d)/2),nrow(d),which(d$gene %in% c('CD3D','CD3E','TRAC','IL7R'))))
  d[chosen,c('cluster','gene','p_val'),drop=FALSE]
}))
cases$cluster<-as.character(cases$cluster);rownames(cases)<-NULL
vectors<-lapply(seq_len(nrow(cases)),function(i) {
  inside<-which(cl$cluster==cases$cluster[[i]])
  outside<-which(cl$cluster!=cases$cluster[[i]])
  list(statistics=as.numeric(X[cases$gene[[i]],c(inside,outside)]),n1=length(inside))
})
rm(obj,X);invisible(gc())
kernel<-function(v)min(2*min(limma::rankSumTestWithCorrelation(
  index=seq_len(v$n1),statistics=v$statistics)),1)
reference<-as.numeric(cases$p_val)
checks<-list();warnings_seen<-character()
for(seed in c(1L,42L,2026L)) {
  set.seed(seed);before<-.Random.seed
  observed<-vapply(vectors,kernel,numeric(1))
  stopifnot(identical(before,.Random.seed),identical(unname(observed),reference))
  checks[[paste0('direct_seed_',seed)]]<-list(seed=seed,workers=1L,
    exact_saved_p_values=TRUE,kernel_did_not_advance_RNG=TRUE)
}
options(future.globals.maxSize=2*1024^3)
for(workers in c(1L,4L))for(seed in c(1L,42L,2026L)) {
  if(workers==1L)plan(sequential) else plan(multicore,workers=workers)
  set.seed(seed)
  observed<-withCallingHandlers(future_sapply(vectors,kernel,future.seed=FALSE),
    warning=function(w)warnings_seen<<-c(warnings_seen,conditionMessage(w)))
  stopifnot(identical(unname(as.numeric(observed)),reference))
  checks[[paste0('future_',workers,'_seed_',seed)]]<-list(seed=seed,workers=workers,exact_saved_p_values=TRUE)
  cat('EXACT_TTU_P_VALUES',workers,seed,nrow(cases),'\n');flush.console()
}
plan(sequential)
proof<-file.path(out,'verification')
write.csv(cases,file.path(proof,'TTU_Wilcoxon_checked_pairs.csv'),row.names=FALSE)
source_text<-c('Seurat WilcoxDETest',deparse(get('WilcoxDETest',asNamespace('Seurat'))),
  'limma rankSumTestWithCorrelation',deparse(limma::rankSumTestWithCorrelation),
  'Seurat FindMarkers.default',deparse(get('FindMarkers.default',asNamespace('Seurat'))))
writeLines(source_text,file.path(proof,'TTU_Wilcoxon_installed_source.txt'))
stopifnot(file.copy(file.path(root,'logs/Rref_score_7351708_4294967294.out'),
  file.path(proof,'TTU_original_SNN_scoring_log.txt'),overwrite=TRUE))
result<-list(status='passed',job=Sys.getenv('SLURM_JOB_ID'),observed_warning_job='7351708',
  cells=ncells,clusters_with_DEG_checked=length(ids),gene_cluster_pairs=nrow(cases),checks=checks,
  input_DEG_sha256=digest(file=file.path(dest,'DEG.rds'),algo='sha256'),
  clusters_sha256=digest(file=file.path(dest,'clusters.csv'),algo='sha256'),
  expression_sha256=digest(file=file.path(prep,'expression_PCA30.rds'),algo='sha256'),
  source_sha256=digest(file=Sys.getenv('DGSCRNA_EXECUTION_SOURCE'),algo='sha256'),
  warnings_during_audit=unique(warnings_seen),
  scope='Deterministic limma rank-sum kernel, all DEG-bearing clusters, first/middle/last saved significant row plus present CD3D/CD3E/TRAC/IL7R rows. This is not a full replay of every tested gene.',
  interpretation='The executed kernel has no stochastic operation; max.cells.per.ident=Inf disables Seurat subsampling. Selected actual saved p-values are exact across three seeds and sequential/four-worker execution. The precise origin of the observed future warning is not established; its original log is retained, and no warning option or scientific parameter was changed.')
write_json(result,file.path(proof,'TTU_wilcox_kernel_audit.json'),auto_unbox=TRUE,pretty=TRUE)
cat('TTU_WILCOX_KERNEL_AUDIT_PASSED',nrow(cases),length(ids),'\n')
