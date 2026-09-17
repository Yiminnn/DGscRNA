#!/usr/bin/env Rscript
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest))
suppressPackageStartupMessages(library(future))
plan(sequential);options(Seurat.object.assay.version='v3')
stopifnot(as.character(packageVersion('dbscan'))=='1.2.6')
i <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
sample <- c('MT-1','MT-2','N-1','N-2','TU-1','TU-2','T-1','T-2')[[i+1L]]
group <- if(i<4L)'MTN' else 'TUT'
base <- '/fs/scratch/PCON0080/yimin/dgscrna/results/hvg_ptc_20260916_v1/ptc_experiments'
dest <- file.path(base,'matched_samples',sample);dir.create(dest,recursive=TRUE,showWarnings=FALSE)
obj <- readRDS(file.path(base,'samples',sample,'RNA_normalized_VST.rds'))
features <- readLines(file.path(base,'prepared',group,'integration_features.txt'))
stopifnot(length(features)==2000L,all(features %in% rownames(obj)))
obj <- ScaleData(obj,features=features,do.center=TRUE,do.scale=TRUE,scale.max=10,verbose=FALSE)
obj <- RunPCA(obj,features=features,npcs=30,seed.use=42,verbose=FALSE)
obj <- RunUMAP(obj,reduction='pca',dims=1:30,n.components=2,n.neighbors=30L,
    umap.method='uwot',metric='cosine',min.dist=0.3,seed.use=42,verbose=FALSE)
write.csv(Embeddings(obj,'pca'),file.path(dest,'PCA30.csv'))
write.csv(Embeddings(obj,'umap'),file.path(dest,'UMAP2.csv'))
writeLines(colnames(obj),file.path(dest,'cells.txt'))
for(space in c('PCA30','UMAP2')) {
  red <- if(space=='PCA30')'pca' else 'umap';z <- Embeddings(obj,red)
  for(method in c('SNN','HDBSCAN_R')) {
    pd <- file.path(dest,paste0(space,'_',method));dir.create(pd,showWarnings=FALSE)
    started <- Sys.time()
    if(method=='SNN') {
      obj <- FindNeighbors(obj,reduction=red,dims=seq_len(ncol(z)),k.param=20,
         compute.SNN=TRUE,prune.SNN=1/15,nn.method='annoy',n.trees=50,annoy.metric='euclidean',
         graph.name=c('ptc_nn','ptc_snn'),verbose=FALSE)
      obj <- FindClusters(obj,graph.name='ptc_snn',resolution=0.5,algorithm=1,
         modularity.fxn=1,n.start=10,n.iter=10,random.seed=0,verbose=FALSE)
      cl <- as.character(Idents(obj));params <- list(resolution=0.5,k=20,algorithm=1,random_seed=0,prune=1/15)
    } else {
      fit <- dbscan::hdbscan(z,minPts=50);cl <- as.character(fit$cluster)
      params <- list(minPts=50,noise=0,implementation='R dbscan1.2.6')
      write.csv(data.frame(cell_id=colnames(obj),membership_prob=fit$membership_prob,
        outlier_score=fit$outlier_scores),file.path(pd,'density_diagnostics.csv'),row.names=FALSE)
    }
    write.csv(data.frame(cell_id=colnames(obj),cluster=cl),file.path(pd,'clusters.csv'),row.names=FALSE)
    write_json(list(sample=sample,group=group,correction='per_sample_NONE',space=space,
      family='matched_R_single',clusterer=method,params=params,n_cells=ncol(obj),
      cluster_sizes=as.list(table(cl)),clusters_sha256=digest(file=file.path(pd,'clusters.csv'),algo='sha256'),
      elapsed_seconds=as.numeric(difftime(Sys.time(),started,units='secs')),annotation_fields_in_fit=FALSE),
      file.path(pd,'cluster_manifest.json'),pretty=TRUE,auto_unbox=TRUE)
    writeLines(digest(file=file.path(pd,'cluster_manifest.json'),algo='sha256'),file.path(pd,'CLUSTER_COMPLETE'))
  }
}
write_json(list(sample=sample,group=group,source_sha256=digest(file=Sys.getenv('PTC_MATCHED_SOURCE'),algo='sha256'),
  feature_sha256=digest(file=file.path(base,'prepared',group,'integration_features.txt'),algo='sha256'),
  n_features=length(features),normalization='same group RNA normalization',
  scaling='Fit separately in this sample; same ScaleData settings',
  seed=42,UMAP=list(neighbors=30,min_dist=0.3,metric='cosine'),
  dbscan='1.2.6',job=Sys.getenv('SLURM_JOB_ID')),
  file.path(dest,'geometry_manifest.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(digest(file=file.path(dest,'geometry_manifest.json'),algo='sha256'),file.path(dest,'GEOMETRY_COMPLETE'))
cat(sample,'matched R single-sample geometry completed\n')
