#!/usr/bin/env Rscript
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest))
suppressPackageStartupMessages(library(future))
plan(sequential)
root<-'/fs/scratch/PCON0080/yimin/dgscrna'
base<-file.path(root,'results/hvg_ptc_20260916_v1/r_reference_campaign_20260917')
args<-commandArgs(trailingOnly=TRUE);group<-args[[1]]
stopifnot(group %in% c('NMT','TTU'))
allprep<-file.path(base,'PTC_ablation',paste0('PTC_',group,'_CCAall'))
stopifnot(file.exists(file.path(allprep,'PREPARED')))
obj<-readRDS(file.path(allprep,'expression_PCA30.rds'))
all_manifest<-fromJSON(file.path(allprep,'prepare_manifest.json'))
all_features<-readLines(file.path(allprep,'DL_features.txt'))
stopifnot(DefaultAssay(obj)=='integrated',identical(all_features,rownames(obj[['integrated']])) )
fixed_features<-readLines(file.path(base,'PTC_ablation',paste0('PTC_',group,'_CCA2000'),'DL_features.txt'))
stopifnot(length(fixed_features)==2000L,all(fixed_features %in% all_features))
fixed_dir<-file.path(base,'PTC_geometry_fixed_inputs',group);dir.create(fixed_dir,recursive=TRUE,showWarnings=FALSE)
fixed_binary<-file.path(fixed_dir,'DL.float32.bin')
con<-file(paste0(fixed_binary,'.part'),'wb')
for(lo in seq.int(1L,ncol(obj),by=256L)) {
  hi<-min(ncol(obj),lo+255L)
  writeBin(as.numeric(as.matrix(obj[['integrated']]@data[fixed_features,lo:hi,drop=FALSE])),con,size=4,endian='little')
}
close(con);stopifnot(file.rename(paste0(fixed_binary,'.part'),fixed_binary))
fixed_hash<-digest(file=fixed_binary,algo='sha256')
writeLines(fixed_features,file.path(fixed_dir,'DL_features.txt'))
fixed_assay<-CreateAssayObject(data=obj[['integrated']]@data[fixed_features,,drop=FALSE])
fixed<-CreateSeuratObject(counts=fixed_assay,assay='fixed_CCA2000',meta.data=obj@meta.data)
for(budget in c('500','1000','2000','3000','5000','all')) {
  unit<-paste0('PTC_',group,'_GEOMETRY',budget,'_FIXED_CCAall_DL2000')
  dest<-file.path(base,'PTC_ablation',unit);dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  if(file.exists(file.path(dest,'PREPARED')))next
  geometry<-if(budget=='all')all_features else readLines(file.path(base,'PTC_ablation',paste0('PTC_',group,'_CCA',budget),'DL_features.txt'))
  stopifnot(all(geometry %in% all_features),length(geometry)>=30L)
  obj<-RunPCA(obj,features=geometry,npcs=30,seed.use=42,verbose=FALSE)
  fixed[['pca']]<-CreateDimReducObject(embeddings=Embeddings(obj,'pca'),key='PC_',assay='fixed_CCA2000')
  saveRDS(fixed,file.path(dest,'expression_PCA30.rds.part'),compress=FALSE)
  stopifnot(file.rename(file.path(dest,'expression_PCA30.rds.part'),file.path(dest,'expression_PCA30.rds')))
  write.csv(data.frame(cell_id=colnames(fixed),batch=fixed$sample_id),file.path(dest,'cells.csv'),row.names=FALSE)
  write.csv(Embeddings(fixed,'pca'),file.path(dest,'PCA30.csv'))
  writeLines(fixed_features,file.path(dest,'DL_features.txt'));writeLines(geometry,file.path(dest,'geometry_features.txt'))
  m<-list(unit=unit,dataset='PTC',group=group,n_cells=ncol(fixed),
    condition=list(group=group,correction='fixed_CCAall_expression',hvg=budget,ablation='geometry_only'),
    correction='fixed_CCAall_expression_geometry_only',assay='fixed_CCA2000',
    features=list(anchor=length(all_features),geometry=length(geometry),scoring=2000L,DL=2000L),
    DL_binary=fixed_binary,DL_binary_sha256=fixed_hash,seed=42L,
    fixed_historical_QC_cells=TRUE,reference_labels_used_for_fitting=FALSE,canonical_cell_order_locked=TRUE,
    fixed_expression_source=file.path(allprep,'prepare_manifest.json'),
    fixed_expression_source_sha256=digest(file=file.path(allprep,'prepare_manifest.json'),algo='sha256'),
    contrast='Only PCA/UMAP geometry features vary; same all-gene CCA fit, scoring matrix, full marker denominators, DL matrix and cell order.',
    limitation='Conditional on an all-shared-gene CCA fit; do not substitute this baseline for the original CCA2000 pipeline.',
    job=Sys.getenv('SLURM_JOB_ID'),source_sha256=digest(file=file.path(root,'handoff/r_reference_campaign_20260917/ptc_geometry_only.R'),algo='sha256'))
  write_json(m,file.path(dest,'prepare_manifest.json'),auto_unbox=TRUE,pretty=TRUE)
  writeLines(capture.output(sessionInfo()),file.path(dest,'sessionInfo.txt'))
  writeLines(digest(file=file.path(dest,'prepare_manifest.json'),algo='sha256'),file.path(dest,'PREPARED'))
  cat('PREPARED',unit,length(geometry),ncol(fixed),'\n');flush.console();gc()
}
