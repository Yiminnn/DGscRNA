stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
extra<-Sys.getenv('DGSCRNA_REFERENCE_R_LIB','')
.libPaths(c(if(nzchar(extra))extra else character(),.Library),include.site=FALSE)
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(jsonlite))
a<-commandArgs(trailingOnly=TRUE);stopifnot(length(a)==2L)
obj<-readRDS(file.path(a[[1]],'expression_PCA30.rds'))
rds<-Embeddings(obj,'umap')
csv<-as.matrix(read.csv(file.path(a[[1]],'UMAP2.csv'),row.names=1,check.names=FALSE,
                       colClasses=c('character','numeric','numeric'),na.strings=NULL))
fr<-dbscan::hdbscan(rds,minPts=50L)
fc<-dbscan::hdbscan(csv,minPts=50L)
stopifnot(identical(fr$cluster,fc$cluster),identical(rownames(rds),rownames(csv)))
for(name in c('rds','csv')) {
  f<-if(name=='rds')fr else fc
  write.csv(data.frame(cell_id=colnames(obj),noise=f$cluster==0L,membership=f$membership_prob),
            file.path(a[[2]],paste0(name,'_membership.csv')),row.names=FALSE)
}
write_json(list(status='same_clusters_different_floating_precision',
  n_cells=nrow(rds),coordinate_max_abs=max(abs(unname(rds)-unname(csv))),
  coordinate_changed=sum(unname(rds)!=unname(csv)),
  membership_max_abs=max(abs(fr$membership_prob-fc$membership_prob)),
  membership_changed=sum(fr$membership_prob!=fc$membership_prob),
  clusters_identical=TRUE,dbscan=as.character(packageVersion('dbscan')),
  Seurat=as.character(packageVersion('Seurat')),R=as.character(getRversion())),
  file.path(a[[2]],'roundtrip.json'),pretty=TRUE,auto_unbox=TRUE,digits=17)
