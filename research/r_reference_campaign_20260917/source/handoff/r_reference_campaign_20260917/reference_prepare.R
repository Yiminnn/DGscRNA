#!/usr/bin/env Rscript
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
# Copy once before evaluating scientific code; queued jobs also enter this guard.
if (!nzchar(Sys.getenv('DGSCRNA_EXECUTION_SOURCE'))) {
  original_script <- '/fs/scratch/PCON0080/yimin/dgscrna/handoff/r_reference_campaign_20260917/reference_prepare.R'
  snapshot_dir <- file.path('/fs/scratch/PCON0080/yimin/dgscrna/results/hvg_ptc_20260916_v1/r_reference_campaign_20260917/execution_sources', Sys.getenv('SLURM_JOB_ID'))
  dir.create(snapshot_dir, recursive=TRUE, showWarnings=FALSE)
  snapshot_script <- file.path(snapshot_dir, 'reference_prepare.R')
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
anchor_workers<-as.integer(Sys.getenv('DGSCRNA_ANCHOR_WORKERS','1'))
stopifnot(anchor_workers>=1L,anchor_workers<=as.integer(Sys.getenv('SLURM_CPUS_PER_TASK','1')))
if(anchor_workers==1L)plan(sequential) else plan(multicore,workers=anchor_workers)
options(Seurat.object.assay.version='v3',future.globals.maxSize=120*1024^3)
root <- '/fs/scratch/PCON0080/yimin/dgscrna'
base <- file.path(root,'results/hvg_ptc_20260916_v1/r_reference_campaign_20260917')
args <- commandArgs(trailingOnly=TRUE)
unit <- if(length(args)) args[[1]] else readLines(file.path(base,'benchmark_units.txt'))[[as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))+1L]]
dataset <- if(length(args)>1L)args[[2]] else if(grepl('^HCL__',unit))'HCL' else unit
inp <- file.path(base,'inputs',dataset,unit)
dest <- file.path(base,'benchmark',unit,'reference_CCA2000');dir.create(dest,recursive=TRUE,showWarnings=FALSE)
im <- fromJSON(file.path(inp,'input_manifest.json'))
if(file.exists(file.path(dest,'PREPARED'))) quit(status=0)
set.seed(42)
log <- function(...) {cat(format(Sys.time()),..., '\n');flush.console();gc()}
save_atomic <- function(x,path) {saveRDS(x,paste0(path,'.part'),compress=FALSE);stopifnot(file.rename(paste0(path,'.part'),path))}
cells <- read.csv(file.path(inp,'cells_fit.csv'),stringsAsFactors=FALSE,check.names=FALSE)
genes <- read.csv(file.path(inp,'genes.csv'),stringsAsFactors=FALSE)$gene
stopifnot(!anyDuplicated(cells$cell_id),!anyDuplicated(genes))
read_bin <- function(name,what,n,size) {con<-file(file.path(inp,name),'rb');on.exit(close(con));readBin(con,what,n=n,size=size,endian='little')}
rna_path <- file.path(dest,'RNA_normalized.rds')
if(file.exists(rna_path)) {
  rna <- readRDS(rna_path)
} else {
  counts <- new('dgCMatrix',x=read_bin('x.bin','numeric',im$nnz,8),i=read_bin('i.bin','integer',im$nnz,4),
       p=read_bin('p.bin','integer',im$n_cells+1L,4),Dim=as.integer(c(im$n_genes,im$n_cells)),Dimnames=list(genes,cells$cell_id))
  rownames(cells) <- cells$cell_id
  rna <- CreateSeuratObject(counts,meta.data=cells,min.cells=0,min.features=0,project=unit)
  rna <- NormalizeData(rna,normalization.method='LogNormalize',scale.factor=10000,verbose=FALSE)
  retained <- Matrix::rowSums(counts>0)>=3L
  rna <- subset(rna,features=rownames(rna)[retained])
  rm(counts);gc()
  save_atomic(rna,rna_path)
}
stopifnot(identical(colnames(rna),cells$cell_id))
batch_sizes <- table(rna$batch)
objects <- SplitObject(rna,split.by='batch')
for(b in names(objects)) {
  objects[[b]] <- FindVariableFeatures(objects[[b]],selection.method='vst',nfeatures=2000,verbose=FALSE)
}
features <- if(length(objects)>1L) SelectIntegrationFeatures(objects,nfeatures=2000) else VariableFeatures(objects[[1]])
stopifnot(length(features)>=30L)
writeLines(features,file.path(dest,'integration_features.txt'))
ndims <- min(30L,min(batch_sizes)-1L,length(features)-1L)
can_integrate <- length(objects)>1L && ndims>=5L
prep_path <- file.path(dest,'expression_PCA30.rds')
if(file.exists(prep_path)) {
  obj <- readRDS(prep_path)
  correction <- readLines(file.path(dest,'correction.txt'))
} else if(can_integrate) {
  log(unit,'CCA anchors',length(objects),'batches',ndims,'dims')
  anchor_path <- file.path(dest,'anchors.rds')
  if(file.exists(anchor_path)) anchors <- readRDS(anchor_path) else {
    anchors <- FindIntegrationAnchors(objects,anchor.features=features,normalization.method='LogNormalize',
      reduction='cca',dims=seq_len(ndims),k.anchor=5,k.filter=min(200L,min(batch_sizes)-1L),
      k.score=min(30L,min(batch_sizes)-1L),max.features=200,nn.method='annoy',n.trees=50,scale=TRUE,verbose=TRUE)
    save_atomic(anchors,anchor_path)
  }
  rm(objects,rna);gc()
  obj <- IntegrateData(anchors,new.assay.name='integrated',normalization.method='LogNormalize',dims=seq_len(ndims),
      k.weight=min(100L,min(batch_sizes)-1L),sd.weight=1,verbose=TRUE)
  rm(anchors);gc()
  DefaultAssay(obj) <- 'integrated';correction <- 'Seurat_CCA_expression'
} else {
  obj <- rna;DefaultAssay(obj)<-'RNA'
  correction <- if(length(objects)==1L)'no_correction_single_batch' else 'CCA_unavailable_batch_too_small_RNA_reference'
  rm(objects,rna);gc()
}
if(!file.exists(prep_path)) {
  VariableFeatures(obj) <- features
  obj <- ScaleData(obj,features=features,do.center=TRUE,do.scale=TRUE,scale.max=10,verbose=FALSE)
  obj <- RunPCA(obj,features=features,npcs=30,seed.use=42,verbose=FALSE)
  save_atomic(obj,prep_path);writeLines(correction,file.path(dest,'correction.txt'))
}
assay <- DefaultAssay(obj)
features <- rownames(obj[[assay]]@data) [rownames(obj[[assay]]@data) %in% features]
writeLines(features,file.path(dest,'DL_features.txt'))
expression <- obj[[assay]]@data[features,,drop=FALSE]
con <- file(file.path(dest,'DL.float32.bin.part'),'wb')
for(lo in seq.int(1L,ncol(expression),by=512L)) {
  hi <- min(ncol(expression),lo+511L)
  writeBin(as.numeric(as.matrix(expression[,lo:hi,drop=FALSE])),con,size=4,endian='little')
}
close(con);stopifnot(file.rename(file.path(dest,'DL.float32.bin.part'),file.path(dest,'DL.float32.bin')))
write.csv(data.frame(cell_id=colnames(obj),batch=obj$batch),file.path(dest,'cells.csv'),row.names=FALSE)
write.csv(Embeddings(obj,'pca'),file.path(dest,'PCA30.csv'))
m <- list(unit=unit,dataset=dataset,n_cells=ncol(obj),input_manifest_sha256=digest(file=file.path(inp,'input_manifest.json'),algo='sha256'),
 input_semantics=im$input_semantics,correction=correction,anchor_dimensions=if(can_integrate)ndims else NULL,
 features=list(anchor=length(features),geometry=length(features),scoring=nrow(obj[[assay]]),DL=length(features)),
 assay=assay,seed=42L,anchor_workers=anchor_workers,DL_binary=file.path(dest,'DL.float32.bin'),DL_binary_sha256=digest(file=file.path(dest,'DL.float32.bin'),algo='sha256'),
 n_batches=length(batch_sizes),job=Sys.getenv('SLURM_JOB_ID'),execution_source=Sys.getenv('DGSCRNA_EXECUTION_SOURCE'),source_sha256=digest(file=Sys.getenv('DGSCRNA_EXECUTION_SOURCE'),algo='sha256'),
 QC='fixed published curated cells; genes detected in >=3 pooled cells; no stochastic re-QC',
 adaptive_parameters='CCA dimension and neighbor caps explicitly reduced only for small batches; see fields and script',
 reference_labels_used_for_fitting=FALSE)
write_json(m,file.path(dest,'prepare_manifest.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(dest,'sessionInfo.txt'))
writeLines(digest(file=file.path(dest,'prepare_manifest.json'),algo='sha256'),file.path(dest,'PREPARED'))
log('PREPARED',unit,correction)
