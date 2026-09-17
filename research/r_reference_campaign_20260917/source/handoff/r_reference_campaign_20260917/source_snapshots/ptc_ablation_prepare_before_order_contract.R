#!/usr/bin/env Rscript
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
# Copy once before evaluating scientific code; queued jobs also enter this guard.
if (!nzchar(Sys.getenv('DGSCRNA_EXECUTION_SOURCE'))) {
  original_script <- '/fs/scratch/PCON0080/yimin/dgscrna/handoff/r_reference_campaign_20260917/ptc_ablation_prepare.R'
  snapshot_dir <- file.path('/fs/scratch/PCON0080/yimin/dgscrna/results/hvg_ptc_20260916_v1/r_reference_campaign_20260917/execution_sources', Sys.getenv('SLURM_JOB_ID'))
  dir.create(snapshot_dir, recursive=TRUE, showWarnings=FALSE)
  snapshot_script <- file.path(snapshot_dir, 'ptc_ablation_prepare.R')
  if (!file.exists(snapshot_script)) stopifnot(file.copy(original_script, snapshot_script))
  Sys.setenv(DGSCRNA_EXECUTION_SOURCE=snapshot_script)
  source(snapshot_script, local=.GlobalEnv)
  quit(status=0)
}

suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest))
suppressPackageStartupMessages(library(future))
plan(sequential);options(Seurat.object.assay.version='v3',future.globals.maxSize=250*1024^3)
root <- '/fs/scratch/PCON0080/yimin/dgscrna'
base <- file.path(root,'results/hvg_ptc_20260916_v1/r_reference_campaign_20260917')
old <- file.path(root,'results/hvg_ptc_20260916_v1/ptc_experiments')
task<-as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
specs<-fromJSON(file.path(root,'handoff/r_reference_campaign_20260917/ptc_ablation_conditions.json'),simplifyVector=FALSE)
spec<-specs[[task+1L]];unit<-spec$unit
dest<-file.path(base,'PTC_ablation',unit);dir.create(dest,recursive=TRUE,showWarnings=FALSE)
if(file.exists(file.path(dest,'PREPARED')))quit(status=0)
set.seed(42)
save_atomic<-function(x,path) {saveRDS(x,paste0(path,'.part'),compress=FALSE);stopifnot(file.rename(paste0(path,'.part'),path))}
sample_groups<-list(NMT=c('MT-1','MT-2','N-1','N-2'),TTU=c('TU-1','TU-2','T-1','T-2'))
group_old<-if(spec$group=='NMT')'MTN' else 'TUT'
if(spec$group=='ALL8') samples<-unlist(sample_groups,use.names=FALSE) else samples<-sample_groups[[spec$group]]
prep_path<-file.path(dest,'expression_PCA30.rds')
if(file.exists(prep_path)) {
  obj<-readRDS(prep_path);features<-readLines(file.path(dest,'DL_features.txt'))
  correction<-readLines(file.path(dest,'correction.txt'))
} else if(spec$correction %in% c('NONE','HARMONY')) {
  # Matched-input geometry contrast: fixed group RNA2000 scoring and DL inputs.
  obj<-readRDS(file.path(old,'prepared',group_old,'NONE_PCA30.rds'))
  features<-readLines(file.path(old,'prepared',group_old,'integration_features.txt'))
  if(spec$correction=='HARMONY') {
    z<-as.matrix(read.csv(file.path(old,'prepared',group_old,'HARMONY_PCA30.csv'),row.names=1,check.names=FALSE))
    stopifnot(setequal(rownames(z),colnames(obj)));z<-z[colnames(obj),,drop=FALSE];colnames(z)<-paste0('PC_',1:30)
    obj[['pca']]<-CreateDimReducObject(embeddings=z,key='PC_',assay='RNA')
  }
  # Restrict scoring to the same RNA features as the CCA2000 expression arm.
  # This isolates expression correction from feature-universe changes.
  obj[['reference_RNA2000']]<-CreateAssayObject(data=obj[['RNA']]@data[features,,drop=FALSE])
  DefaultAssay(obj)<-'reference_RNA2000'
  correction<-paste0(spec$correction,'_RNA2000_scoring_DL')
} else {
  objects<-list()
  for(sample in samples) {
    o<-readRDS(file.path(old,'samples',sample,'RNA_normalized_VST.rds'))
    eligible<-readLines(file.path(old,'samples',sample,'geometry_genes.txt'))
    o<-subset(o,features=eligible)
    # Remove archived annotations from the fitting object.
    allowed<-intersect(c('orig.ident','nCount_RNA','nFeature_RNA','sample_id','percent.mt'),colnames(o@meta.data))
    o@meta.data<-o@meta.data[,allowed,drop=FALSE]
    o$sample_id<-sample
    o<-FindVariableFeatures(o,selection.method='vst',nfeatures=if(spec$hvg=='all')nrow(o) else as.integer(spec$hvg),verbose=FALSE)
    objects[[sample]]<-o;rm(o);gc()
  }
  if(spec$hvg=='all') {
    features<-Reduce(intersect,lapply(objects,rownames))
  } else features<-SelectIntegrationFeatures(objects,nfeatures=as.integer(spec$hvg))
  stopifnot(length(features)>=30L)
  writeLines(features,file.path(dest,'DL_features.txt'))
  anchorpath<-file.path(dest,'anchors.rds')
  if(file.exists(anchorpath)) anchors<-readRDS(anchorpath) else {
    anchors<-FindIntegrationAnchors(objects,anchor.features=features,normalization.method='LogNormalize',reduction='cca',
       dims=1:30,k.anchor=5,k.filter=200,k.score=30,max.features=200,nn.method='annoy',n.trees=50,scale=TRUE,verbose=TRUE)
    save_atomic(anchors,anchorpath)
  }
  rm(objects);gc()
  obj<-IntegrateData(anchors,new.assay.name='integrated',normalization.method='LogNormalize',dims=1:30,
        k.weight=100,sd.weight=1,verbose=TRUE)
  rm(anchors);gc();DefaultAssay(obj)<-'integrated'
  obj<-ScaleData(obj,features=features,do.center=TRUE,do.scale=TRUE,scale.max=10,verbose=FALSE)
  obj<-RunPCA(obj,features=features,npcs=30,seed.use=42,verbose=FALSE)
  correction<-paste0('fresh_Seurat_CCA_',spec$group,'_',spec$hvg)
}
assay<-DefaultAssay(obj)
features<-rownames(obj[[assay]]@data)[rownames(obj[[assay]]@data) %in% features]
writeLines(features,file.path(dest,'DL_features.txt'));writeLines(correction,file.path(dest,'correction.txt'))
if(!file.exists(prep_path))save_atomic(obj,prep_path)
con<-file(file.path(dest,'DL.float32.bin.part'),'wb')
for(lo in seq.int(1L,ncol(obj),by=256L)) {
  hi<-min(ncol(obj),lo+255L)
  writeBin(as.numeric(as.matrix(obj[[assay]]@data[features,lo:hi,drop=FALSE])),con,size=4,endian='little')
}
close(con);stopifnot(file.rename(file.path(dest,'DL.float32.bin.part'),file.path(dest,'DL.float32.bin')))
write.csv(data.frame(cell_id=colnames(obj),batch=obj$sample_id),file.path(dest,'cells.csv'),row.names=FALSE)
write.csv(Embeddings(obj,'pca'),file.path(dest,'PCA30.csv'))
m<-list(unit=unit,dataset='PTC',group=spec$group,condition=spec,n_cells=ncol(obj),correction=correction,
 features=list(anchor=if(spec$correction=='CCA')length(features) else 0L,geometry=length(features),scoring=nrow(obj[[assay]]),DL=length(features)),
 assay=assay,seed=42L,DL_binary=file.path(dest,'DL.float32.bin'),DL_binary_sha256=digest(file=file.path(dest,'DL.float32.bin'),algo='sha256'),
 reference_labels_used_for_fitting=FALSE,fixed_historical_QC_cells=TRUE,
 no_HVG_definition=if(spec$hvg=='all')'all genes detected in >=3 cells in every included sample; no variance ranking' else NULL,
 correction_comparison='CCA versus RNA-based correction changes expression as well as geometry; NONE versus HARMONY holds score/DL RNA2000 fixed',
 job=Sys.getenv('SLURM_JOB_ID'),execution_source=Sys.getenv('DGSCRNA_EXECUTION_SOURCE'),source_sha256=digest(file=Sys.getenv('DGSCRNA_EXECUTION_SOURCE'),algo='sha256'))
write_json(m,file.path(dest,'prepare_manifest.json'),auto_unbox=TRUE,pretty=TRUE)
writeLines(capture.output(sessionInfo()),file.path(dest,'sessionInfo.txt'))
writeLines(digest(file=file.path(dest,'prepare_manifest.json'),algo='sha256'),file.path(dest,'PREPARED'))
cat('PREPARED',unit,ncol(obj),length(features),'\n')
