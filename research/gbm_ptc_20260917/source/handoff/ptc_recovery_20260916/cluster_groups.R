#!/usr/bin/env Rscript
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest))
suppressPackageStartupMessages(library(future))
plan(sequential);options(Seurat.object.assay.version='v3')
args <- commandArgs(trailingOnly=TRUE)
task <- if(length(args)) as.integer(args[[1]]) else as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
group <- c('MTN','TUT')[[task %/% 3L+1L]]
correction <- c('NONE','CCA','HARMONY')[[task %% 3L+1L]]
root <- '/fs/scratch/PCON0080/yimin/dgscrna'
base <- file.path(root,'results/hvg_ptc_20260916_v1/ptc_experiments')
prep <- file.path(base,'prepared',group)
dest <- file.path(base,'pooled',group,correction)
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
stopifnot(file.exists(file.path(prep,'PREPARATION_COMPLETE')))
if(correction=='HARMONY') stopifnot(file.exists(file.path(prep,'HARMONY_COMPLETE')))
embpath <- file.path(prep,paste0(correction,'_PCA30.csv'))
emb <- as.matrix(read.csv(embpath,row.names=1,check.names=FALSE))
rna <- readRDS(file.path(prep,'RNA_normalized.rds'))
stopifnot(setequal(rownames(emb),colnames(rna)))
emb <- emb[colnames(rna),,drop=FALSE]
colnames(emb) <- paste0('PTCPC_',1:30)
rna[['ptcpc']] <- CreateDimReducObject(embeddings=emb,key='PTCPC_',assay='RNA')
umap_path <- file.path(dest,'UMAP2.csv')
if(file.exists(umap_path)) {
  umap <- as.matrix(read.csv(umap_path,row.names=1,check.names=FALSE))[colnames(rna),,drop=FALSE]
  colnames(umap) <- paste0('PTCUMAP_',1:2)
  rna[['ptcumap']] <- CreateDimReducObject(embeddings=umap,key='PTCUMAP_',assay='RNA')
} else {
  rna <- RunUMAP(rna,reduction='ptcpc',dims=1:30,n.components=2,n.neighbors=30L,
     umap.method='uwot',metric='cosine',min.dist=0.3,seed.use=42,
     reduction.name='ptcumap',reduction.key='PTCUMAP_',verbose=TRUE)
  write.csv(Embeddings(rna,'ptcumap'),umap_path)
}
write.csv(emb,file.path(dest,'PCA30.csv'))
records <- list()
for(space in c('PCA30','UMAP2')) {
  red <- if(space=='PCA30') 'ptcpc' else 'ptcumap'
  z <- Embeddings(rna,red)
  for(method in c('SNN','HDBSCAN_R')) {
    pd <- file.path(dest,paste0(space,'_',method));dir.create(pd,showWarnings=FALSE)
    if(file.exists(file.path(pd,'CLUSTER_COMPLETE'))) next
    cat(format(Sys.time(),tz='UTC'),group,correction,space,method,'\n');flush.console()
    started <- Sys.time()
    if(method=='SNN') {
      rna <- FindNeighbors(rna,reduction=red,dims=seq_len(ncol(z)),k.param=20,
          compute.SNN=TRUE,prune.SNN=1/15,nn.method='annoy',n.trees=50,annoy.metric='euclidean',
          graph.name=c('ptc_nn','ptc_snn'),verbose=FALSE)
      rna <- FindClusters(rna,graph.name='ptc_snn',resolution=0.5,algorithm=1,
          modularity.fxn=1,n.start=10,n.iter=10,random.seed=0,verbose=FALSE)
      cl <- as.character(Idents(rna))
      params <- list(resolution=0.5,algorithm=1,k=20,random_seed=0,prune=1/15)
    } else {
      fit <- dbscan::hdbscan(z,minPts=50)
      cl <- as.character(fit$cluster)
      write.csv(data.frame(cell_id=rownames(z),membership_prob=fit$membership_prob,
                          outlier_score=fit$outlier_scores),file.path(pd,'density_diagnostics.csv'),row.names=FALSE)
      params <- list(minPts=50,noise=0,implementation=paste0('R dbscan ',packageVersion('dbscan')))
    }
    write.csv(data.frame(cell_id=colnames(rna),cluster=cl),file.path(pd,'clusters.csv'),row.names=FALSE)
    m <- list(group=group,correction=correction,space=space,clusterer=method,params=params,
      job=Sys.getenv('SLURM_JOB_ID'),n_cells=length(cl),cluster_sizes=as.list(table(cl)),
      input_embedding_sha256=digest(file=embpath,algo='sha256'),
      clusters_sha256=digest(file=file.path(pd,'clusters.csv'),algo='sha256'),
      elapsed_seconds=as.numeric(difftime(Sys.time(),started,units='secs')),
      annotation_fields_in_fit=FALSE,noise_scored_as_observed_cluster=TRUE)
    write_json(m,file.path(pd,'cluster_manifest.json'),pretty=TRUE,auto_unbox=TRUE)
    writeLines(digest(file=file.path(pd,'cluster_manifest.json'),algo='sha256'),file.path(pd,'CLUSTER_COMPLETE'))
    records[[paste(space,method)]] <- m
  }
}
write_json(list(group=group,correction=correction,job=Sys.getenv('SLURM_JOB_ID'),
  source_sha256=digest(file=Sys.getenv('PTC_CLUSTER_SOURCE'),algo='sha256'),
  UMAP=list(neighbors=30,min_dist=0.3,metric='cosine',seed=42),records=records),
  file.path(dest,'geometry_manifest.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(digest(file=file.path(dest,'geometry_manifest.json'),algo='sha256'),file.path(dest,'GEOMETRY_COMPLETE'))
