#!/usr/bin/env Rscript
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest))
root<-'/fs/scratch/PCON0080/yimin/dgscrna'
.libPaths(c(file.path(root,'results/hvg_ptc_20260916_v1/ptc_experiments/vendor_R_dbscan64'),.libPaths()))
args<-commandArgs(trailingOnly=TRUE);mode<-args[[1]];cfg<-fromJSON(args[[2]])
gm<-fromJSON(file.path(cfg$geometry,'manifest.json'))
hash<-function(p)digest(file=p,algo='sha256')
stopifnot(readLines(file.path(cfg$geometry,'COMPLETE'))==hash(file.path(cfg$geometry,'manifest.json')))
cells<-read.csv(file.path(cfg$geometry,'cells.csv'),stringsAsFactors=FALSE)$cell_id
readmat<-function(p,n,d) {con<-file(p,'rb');on.exit(close(con));matrix(readBin(con,'numeric',n=n*d,size=8,endian='little'),nrow=n,byrow=TRUE)}
if(mode=='umap') {
  stopifnot(hash(gm$binary)==gm$binary_sha256)
  x<-readmat(gm$binary,gm$n_cells,gm$n_features)
  set.seed(42)
  z<-uwot::umap(x,n_neighbors=30L,n_components=2L,metric='cosine',min_dist=0.3,
    spread=1,learning_rate=1,repulsion_strength=1,negative_sample_rate=5,
    init='spectral',n_threads=1,n_sgd_threads=1,verbose=TRUE)
  stopifnot(all(is.finite(z)),nrow(z)==length(cells))
  write.csv(data.frame(cell_id=cells,x=z[,1],y=z[,2]),cfg$embedding,row.names=FALSE)
  write_json(list(implementation='uwot::umap',version=as.character(packageVersion('uwot')),
    seed=42,n_neighbors=30,n_components=2,metric='cosine',min_dist=0.3,
    input='direct native scaled-HVG; no PCA',threads=1),paste0(cfg$embedding,'.params.json'),pretty=TRUE,auto_unbox=TRUE)
} else if(mode=='hdbscan') {
  if(cfg$space=='noDR') {
    stopifnot(hash(gm$binary)==gm$binary_sha256)
    z<-readmat(gm$binary,gm$n_cells,gm$n_features)
  } else {
    e<-read.csv(cfg$embedding,stringsAsFactors=FALSE)
    stopifnot(identical(e$cell_id,cells));z<-as.matrix(e[,c('x','y')])
  }
  fit<-dbscan::hdbscan(z,minPts=50L)
  dir.create(cfg$hdbscan_dest,recursive=TRUE,showWarnings=FALSE)
  write.csv(data.frame(cell_id=cells,cluster=fit$cluster),file.path(cfg$hdbscan_dest,'clusters.csv'),row.names=FALSE)
  write.csv(data.frame(cell_id=cells,noise=fit$cluster==0L,membership=fit$membership_prob),file.path(cfg$hdbscan_dest,'density_diagnostics.csv'),row.names=FALSE)
  write_json(list(implementation='dbscan::hdbscan',version=as.character(packageVersion('dbscan')),
    minPts=50,noise_label=0,noise_scored_as_cluster=TRUE),file.path(cfg$hdbscan_dest,'cluster_params.json'),pretty=TRUE,auto_unbox=TRUE)
} else stop('Unknown mode')
