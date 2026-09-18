#!/usr/bin/env Rscript
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[[1]])
Sys.setenv(DGSCRNA_EXECUTION_SOURCE=script)
root <- '/fs/scratch/PCON0080/yimin/dgscrna'
base <- file.path(root,'results/hvg_ptc_20260916_v1/paper_claim_validation_20260917')
.libPaths(c(file.path(root,'results/hvg_ptc_20260916_v1/ptc_experiments/vendor_R_dbscan64'),.libPaths()))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest))
suppressPackageStartupMessages(library(future))
deg_workers<-min(4L,as.integer(Sys.getenv('SLURM_CPUS_PER_TASK','4')))
plan(multicore,workers=deg_workers)
options(future.globals.maxSize=32*1024^3)
args <- commandArgs(trailingOnly=TRUE)
unit <- args[[1]]
only_route<-Sys.getenv('DGSCRNA_ONLY_ROUTE','')
all_routes<-c('PCA30_SNN','PCA30_HDBSCAN_R','UMAP2_SNN','UMAP2_HDBSCAN_R')
stopifnot(!nzchar(only_route) || only_route %in% all_routes)
prep <- if(length(args)>1L) args[[2]] else file.path(base,'benchmark',unit,'reference_CCA2000')
stopifnot(file.exists(file.path(prep,'PREPARED')))
stopifnot(readLines(file.path(prep,'PREPARED'))==digest(file=file.path(prep,'prepare_manifest.json'),algo='sha256'))
pm <- fromJSON(file.path(prep,'prepare_manifest.json'))
# Lower fork concurrency for the all-gene assay; same DEG statistics.
if(pm$features$scoring>10000L) {
  deg_workers<-min(as.integer(Sys.getenv('DGSCRNA_DEG_WORKERS','2')),as.integer(Sys.getenv('SLURM_CPUS_PER_TASK','2')))
  stopifnot(deg_workers>=1L,deg_workers<=4L)
  plan(multicore,workers=deg_workers)
}
stopifnot(digest(file=file.path(prep,'expression_PCA30.rds'),algo='sha256')==pm$expression_sha256)
obj <- readRDS(file.path(prep,'expression_PCA30.rds'))
cells <- read.csv(file.path(prep,'cells.csv'),stringsAsFactors=FALSE)
stopifnot(identical(colnames(obj),cells$cell_id))
publish <- function(path,writer) {
  tmp<-paste0(path,'.part.',Sys.getenv('SLURM_JOB_ID'))
  writer(tmp)
  stopifnot(file.rename(tmp,path))
}
csv_atomic <- function(x,path,...) publish(path,function(p)write.csv(x,p,...))
rds_atomic <- function(x,path) publish(path,function(p)saveRDS(x,p))
read_cache <- function(path,reader,valid) {
  if(!file.exists(path))return(NULL)
  tryCatch({x<-reader(path);stopifnot(isTRUE(valid(x)));x},error=function(e) {
    preserved<-tempfile(paste0(basename(path),'.invalid_'),tmpdir=dirname(path))
    stopifnot(file.rename(path,preserved))
    event<-list(path=path,preserved_as=preserved,error=conditionMessage(e),job=Sys.getenv('SLURM_JOB_ID'))
    cat(toJSON(event,auto_unbox=TRUE),'\n',file=file.path(prep,paste0('cache_recovery_',Sys.getenv('SLURM_JOB_ID'),'.jsonl')),append=TRUE)
    cat('INVALID_CACHE_PRESERVED',path,preserved,conditionMessage(e),'\n');flush.console()
    NULL
  })
}
marker_path <- file.path(base,'markers','libraries.json')
libraries <- fromJSON(marker_path,simplifyVector=FALSE)
libraries <- lapply(libraries,function(l)lapply(l,unlist,use.names=FALSE))
stopifnot(length(libraries)>1L)
umap_path <- file.path(prep,'UMAP2.csv')
embedding <- read_cache(umap_path,function(p)as.matrix(read.csv(p,row.names=1,check.names=FALSE)),
  function(x)identical(rownames(x),colnames(obj)) && ncol(x)==2L && all(is.finite(x)))
if(!is.null(embedding)) {
  colnames(embedding)<-paste0('UMAP_',1:2)
  obj[['umap']] <- CreateDimReducObject(embeddings=embedding,key='UMAP_',assay=DefaultAssay(obj))
} else {
  obj <- RunUMAP(obj,reduction='pca',dims=1:30,n.components=2,n.neighbors=30L,
       umap.method='uwot',metric='cosine',min.dist=0.3,seed.use=42,verbose=FALSE)
  csv_atomic(Embeddings(obj,'umap'),umap_path)
}
assay <- DefaultAssay(obj)
legacy_mean <- function(x) log2(Matrix::rowMeans(expm1(x))+1)
for(space in c('PCA30','UMAP2')) for(method in c('SNN','HDBSCAN_R')) {
  route <- paste(space,method,sep='_')
  if(nzchar(only_route) && route!=only_route)next
  dest <- file.path(prep,route);dir.create(dest,showWarnings=FALSE)
  if(file.exists(file.path(dest,'SCORE_COMPLETE'))) next
  started <- Sys.time()
  red <- if(space=='PCA30')'pca' else 'umap'
  z <- Embeddings(obj,red)
  cf <- read_cache(file.path(dest,'clusters.csv'),function(p)read.csv(p,stringsAsFactors=FALSE),
    function(x)all(c('cell_id','cluster') %in% names(x)) && identical(x$cell_id,colnames(obj)) && !anyNA(x$cluster))
  if(!is.null(cf)) {
    cl<-as.character(cf$cluster)
  } else {
    if(method=='SNN') {
      obj <- FindNeighbors(obj,reduction=red,dims=seq_len(ncol(z)),k.param=20L,compute.SNN=TRUE,
        prune.SNN=1/15,nn.method='annoy',n.trees=50,annoy.metric='euclidean',graph.name=c('ref_nn','ref_snn'),verbose=FALSE)
      obj <- FindClusters(obj,graph.name='ref_snn',resolution=0.5,algorithm=1,modularity.fxn=1,
        n.start=10,n.iter=10,random.seed=0,verbose=FALSE)
      cl <- as.character(Idents(obj))
    } else {
      fit <- dbscan::hdbscan(z,minPts=50L)
      cl <- as.character(fit$cluster)
      csv_atomic(data.frame(cell_id=colnames(obj),noise=fit$cluster==0L,membership=fit$membership_prob),
        file.path(dest,'density_diagnostics.csv'),row.names=FALSE)
    }
    csv_atomic(data.frame(cell_id=colnames(obj),cluster=cl),file.path(dest,'clusters.csv'),row.names=FALSE)
  }
  ids <- sort(unique(cl));Idents(obj)<-factor(cl,levels=ids)
  cat(format(Sys.time()),unit,route,'clusters',length(ids),'DEG\n');flush.console()
  degpath <- file.path(dest,'DEG.rds')
  deg<-read_cache(degpath,readRDS,function(x)is.data.frame(x) && all(c('cluster','gene','avg_log2FC') %in% names(x)) &&
    all(as.character(x$cluster) %in% ids) && !anyNA(x$gene) && all(is.finite(x$avg_log2FC)))
  if(is.null(deg)) {
    # Release completed clustering temporaries before forked DEG workers inherit the heap.
    invisible(gc())
    if(length(ids)>1L) {
      deg<-FindAllMarkers(obj,assay=assay,slot='data',test.use='wilcox_limma',
        logfc.threshold=0.25,min.pct=0.1,min.diff.pct=-Inf,only.pos=FALSE,
        max.cells.per.ident=Inf,random.seed=1,min.cells.feature=3,min.cells.group=3,
        pseudocount.use=1,mean.fxn=legacy_mean,fc.name='avg_log2FC',base=2,
        return.thresh=0.01,densify=FALSE,verbose=FALSE)
    } else deg <- data.frame()
    if(nrow(deg)==0L) deg<-data.frame(cluster=character(),gene=character(),avg_log2FC=numeric())
    stopifnot(all(is.finite(deg$avg_log2FC)))
    rds_atomic(deg,degpath)
  }
  # The arithmetic is the archived density score. Explicit observed-cluster
  # indexing prevents the legacy 0..K-1 assumption from losing a noiseless HDBSCAN group.
  initial<-data.frame(cell_id=colnames(obj),stringsAsFactors=FALSE)
  arms<-list();callrows<-list();retention<-list();score_matrices<-list()
  for(i in seq_along(libraries)) {
    lib<-names(libraries)[[i]];panels<-libraries[[i]]
    S<-matrix(0,length(panels),length(ids),dimnames=list(names(panels),ids))
    for(id in ids) {
      ds<-deg[as.character(deg$cluster)==id & deg$avg_log2FC>1,,drop=FALSE]
      fc<-setNames(ds$avg_log2FC,ds$gene)
      for(p in names(panels)) {
        # Preserve the original DEG-row summation order as well as its formula.
        S[p,id]<-sum(fc[names(fc) %in% panels[[p]]])/length(panels[[p]])
        if(length(panels[[p]])<=1L)S[p,id]<-S[p,id]*0.8
      }
    }
    score_matrices[[lib]]<-S;mx<-apply(S,2,max)
    retention[[lib]]<-data.frame(library=lib,panel=names(panels),denominator=lengths(panels),
      retained=vapply(panels,function(g)length(intersect(g,rownames(obj[[assay]]))),integer(1)))
    for(cut in c('none','mean','0.5')) {
      calls<-vapply(ids,function(id) {
        winners<-rownames(S)[S[,id]==mx[[id]]]
        answer<-if(length(winners)==1L) winners else 'Undecided'
        threshold<-if(cut=='mean') mean(mx) else if(cut=='0.5')0.5 else -Inf
        if(mx[[id]]<threshold)answer<-'Undecided'
        answer
      },character(1))
      aid<-sprintf('L%02d_%s',i-1L,if(cut=='0.5')'p050' else cut)
      initial[[aid]]<-unname(calls[cl]);stopifnot(!anyNA(initial[[aid]]))
      arms[[aid]]<-list(arm_id=aid,library=lib,cutoff=cut,seed_column=aid)
      callrows[[aid]]<-data.frame(arm_id=aid,library=lib,cutoff=cut,cluster=ids,initial=unname(calls),max_score=as.numeric(mx))
    }
  }
  write_gz<-function(x,path) publish(path,function(p) {con<-gzfile(p,'wt');write.csv(x,con,row.names=FALSE);close(con)})
  write_gz(initial,file.path(dest,'initial_calls.csv.gz'))
  write_gz(do.call(rbind,callrows),file.path(dest,'cluster_calls.csv.gz'))
  write_gz(do.call(rbind,retention),file.path(dest,'marker_retention.csv.gz'))
  rds_atomic(score_matrices,file.path(dest,'density_scores.rds'))
  csv_atomic(cells,file.path(dest,'cells.csv'),row.names=FALSE)
  m<-list(status='score_complete_DL_pending',dataset=pm$dataset,unit=unit,route=route,budget=pm$budget,seed=pm$seed,
    geometry=pm$correction,assay=assay,n_cells=ncol(obj),n_clusters=length(ids),
    scoring_features=nrow(obj[[assay]]),DL_features=pm$features$DL,DL_binary=pm$DL_binary,
    DL_binary_sha256=pm$DL_binary_sha256,DL_input_description=paste(pm$correction,assay,pm$features$DL,'features'),
    marker_source_sha256=digest(file=marker_path,algo='sha256'),
    initial_sha256=digest(file=file.path(dest,'initial_calls.csv.gz'),algo='sha256'),
    execution_source=Sys.getenv('DGSCRNA_EXECUTION_SOURCE'),source_sha256=digest(file=Sys.getenv('DGSCRNA_EXECUTION_SOURCE'),algo='sha256'),
    original_density_arithmetic=TRUE,cluster_indexing='observed labels; equivalent to legacy contiguous zero-based indexing',
    noise_zero_scored_as_cluster=TRUE,dbscan_version=as.character(packageVersion('dbscan')),arms=arms,
    reference_labels_used_for_fit=FALSE,DEG_workers=deg_workers,job=Sys.getenv('SLURM_JOB_ID'),
    elapsed_seconds=as.numeric(difftime(Sys.time(),started,units='secs')))
  publish(file.path(dest,'score_manifest.json'),function(p)write_json(m,p,auto_unbox=TRUE,pretty=TRUE))
  publish(file.path(dest,'SCORE_COMPLETE'),function(p)writeLines(digest(file=file.path(dest,'score_manifest.json'),algo='sha256'),p))
  cat('SCORE_COMPLETE',unit,route,length(arms),'arms\n');flush.console()
}
if(all(file.exists(file.path(prep,all_routes,'SCORE_COMPLETE'))))
  publish(file.path(prep,'GRID_SCORED'),function(p)writeLines('four original clustering branches scored',p))
