#!/usr/bin/env Rscript
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
root <- '/fs/scratch/PCON0080/yimin/dgscrna'
base <- file.path(root,'results/hvg_ptc_20260916_v1/r_reference_campaign_20260917')
.libPaths(c(file.path(root,'results/hvg_ptc_20260916_v1/ptc_experiments/vendor_R_dbscan64'),.libPaths()))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest))
suppressPackageStartupMessages(library(future))
plan(multicore,workers=as.integer(Sys.getenv('SLURM_CPUS_PER_TASK','4')))
options(future.globals.maxSize=32*1024^3)
args <- commandArgs(trailingOnly=TRUE)
unit <- args[[1]]
prep <- if(length(args)>1L) args[[2]] else file.path(base,'benchmark',unit,'reference_CCA2000')
stopifnot(file.exists(file.path(prep,'PREPARED')))
pm <- fromJSON(file.path(prep,'prepare_manifest.json'))
obj <- readRDS(file.path(prep,'expression_PCA30.rds'))
cells <- read.csv(file.path(prep,'cells.csv'),stringsAsFactors=FALSE)
stopifnot(identical(colnames(obj),cells$cell_id))
marker_path <- file.path(base,'markers',if(startsWith(unit,'PTC_'))'PTC_original17.json' else paste0(unit,'.json'))
libraries <- fromJSON(marker_path,simplifyVector=FALSE)
libraries <- lapply(libraries,function(l)lapply(l,unlist,use.names=FALSE))
stopifnot(length(libraries)>1L)
umap_path <- file.path(prep,'UMAP2.csv')
if(file.exists(umap_path)) {
  embedding <- as.matrix(read.csv(umap_path,row.names=1,check.names=FALSE))
  stopifnot(identical(rownames(embedding),colnames(obj)))
  colnames(embedding)<-paste0('UMAP_',1:2)
  obj[['umap']] <- CreateDimReducObject(embeddings=embedding,key='UMAP_',assay=DefaultAssay(obj))
} else {
  obj <- RunUMAP(obj,reduction='pca',dims=1:30,n.components=2,n.neighbors=30L,
       umap.method='uwot',metric='cosine',min.dist=0.3,seed.use=42,verbose=FALSE)
  write.csv(Embeddings(obj,'umap'),umap_path)
}
assay <- DefaultAssay(obj)
legacy_mean <- function(x) log2(Matrix::rowMeans(expm1(x))+1)
for(space in c('PCA30','UMAP2')) for(method in c('SNN','HDBSCAN_R')) {
  route <- paste(space,method,sep='_')
  dest <- file.path(prep,route);dir.create(dest,showWarnings=FALSE)
  if(file.exists(file.path(dest,'SCORE_COMPLETE'))) next
  started <- Sys.time()
  red <- if(space=='PCA30')'pca' else 'umap'
  z <- Embeddings(obj,red)
  if(file.exists(file.path(dest,'clusters.csv'))) {
    cf <- read.csv(file.path(dest,'clusters.csv'),stringsAsFactors=FALSE)
    stopifnot(identical(cf$cell_id,colnames(obj)));cl<-as.character(cf$cluster)
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
      write.csv(data.frame(cell_id=colnames(obj),noise=fit$cluster==0L,membership=fit$membership_prob),
        file.path(dest,'density_diagnostics.csv'),row.names=FALSE)
    }
    write.csv(data.frame(cell_id=colnames(obj),cluster=cl),file.path(dest,'clusters.csv'),row.names=FALSE)
  }
  ids <- sort(unique(cl));Idents(obj)<-factor(cl,levels=ids)
  cat(format(Sys.time()),unit,route,'clusters',length(ids),'DEG\n');flush.console()
  degpath <- file.path(dest,'DEG.rds')
  if(file.exists(degpath)) deg<-readRDS(degpath) else {
    if(length(ids)>1L) {
      deg<-FindAllMarkers(obj,assay=assay,slot='data',test.use='wilcox_limma',
        logfc.threshold=0.25,min.pct=0.1,min.diff.pct=-Inf,only.pos=FALSE,
        max.cells.per.ident=Inf,random.seed=1,min.cells.feature=3,min.cells.group=3,
        pseudocount.use=1,mean.fxn=legacy_mean,fc.name='avg_log2FC',base=2,
        return.thresh=0.01,densify=FALSE,verbose=FALSE)
    } else deg <- data.frame()
    if(nrow(deg)==0L) deg<-data.frame(cluster=character(),gene=character(),avg_log2FC=numeric())
    stopifnot(all(is.finite(deg$avg_log2FC)))
    saveRDS(deg,degpath)
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
        hits<-intersect(panels[[p]],names(fc))
        S[p,id]<-sum(fc[hits])/length(panels[[p]])
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
  write_gz<-function(x,path) {con<-gzfile(path,'wt');write.csv(x,con,row.names=FALSE);close(con)}
  write_gz(initial,file.path(dest,'initial_calls.csv.gz'))
  write_gz(do.call(rbind,callrows),file.path(dest,'cluster_calls.csv.gz'))
  write_gz(do.call(rbind,retention),file.path(dest,'marker_retention.csv.gz'))
  saveRDS(score_matrices,file.path(dest,'density_scores.rds'))
  write.csv(cells,file.path(dest,'cells.csv'),row.names=FALSE)
  m<-list(status='score_complete_DL_pending',dataset=pm$dataset,unit=unit,route=route,
    geometry=pm$correction,assay=assay,n_cells=ncol(obj),n_clusters=length(ids),
    scoring_features=nrow(obj[[assay]]),DL_features=pm$features$DL,DL_binary=pm$DL_binary,
    DL_binary_sha256=pm$DL_binary_sha256,DL_input_description=paste(pm$correction,assay,pm$features$DL,'features'),
    marker_source_sha256=digest(file=marker_path,algo='sha256'),
    initial_sha256=digest(file=file.path(dest,'initial_calls.csv.gz'),algo='sha256'),
    source_sha256=digest(file=file.path(root,'handoff/r_reference_campaign_20260917/reference_score.R'),algo='sha256'),
    original_density_arithmetic=TRUE,cluster_indexing='observed labels; equivalent to legacy contiguous zero-based indexing',
    noise_zero_scored_as_cluster=TRUE,dbscan_version=as.character(packageVersion('dbscan')),arms=arms,
    reference_labels_used_for_fit=FALSE,job=Sys.getenv('SLURM_JOB_ID'),
    elapsed_seconds=as.numeric(difftime(Sys.time(),started,units='secs')))
  write_json(m,file.path(dest,'score_manifest.json'),auto_unbox=TRUE,pretty=TRUE)
  writeLines(digest(file=file.path(dest,'score_manifest.json'),algo='sha256'),file.path(dest,'SCORE_COMPLETE'))
  cat('SCORE_COMPLETE',unit,route,length(arms),'arms\n');flush.console()
}
writeLines('four original clustering branches scored',file.path(prep,'GRID_SCORED'))
