#!/usr/bin/env Rscript
# Fixed-QC R preprocessing and original-style CCA in the two user-specified groups.
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
args <- commandArgs(trailingOnly=TRUE)
group <- if(length(args)) args[[1]] else c('MTN','TUT')[[as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))+1L]]
stopifnot(group %in% c('MTN','TUT'))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest))
suppressPackageStartupMessages(library(future))
plan(sequential)
options(Seurat.object.assay.version='v3',future.globals.maxSize=150*1024^3)
set.seed(42)
root <- '/fs/scratch/PCON0080/yimin/dgscrna'
recovery <- file.path(root,'results/hvg_ptc_20260916_v1/ptc_recovery')
base <- file.path(root,'results/hvg_ptc_20260916_v1/ptc_experiments')
dest <- file.path(base,'prepared',group)
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
stopifnot(file.exists(file.path(recovery,'R_baseline_replay/COMPLETE')))
source_path <- file.path(recovery,'fixed_QC_inputs',paste0(group,'.counts.rds'))
protocol_path <- file.path(root,'handoff/ptc_recovery_20260916/CONTROLLED_PROTOCOL.md')
input_sha <- digest(file=source_path,algo='sha256')
protocol_sha <- digest(file=protocol_path,algo='sha256')
start <- Sys.time()
log_step <- function(text) {
  cat(format(Sys.time(),tz='UTC'),text,'\n');flush.console()
  memory <- readLines('/proc/self/status')
  cat(paste(memory[grepl('VmRSS|VmHWM',memory)],collapse='; '),'\n');flush.console()
}
save_checkpoint <- function(obj,path) {
  tmp <- paste0(path,'.part');saveRDS(obj,tmp,compress=FALSE)
  stopifnot(file.rename(tmp,path))
}
write_binary <- function(mat,path) {
  # R feature×cell column-major order is Python cell×feature C order.
  con <- file(paste0(path,'.part'),'wb')
  on.exit(close(con))
  chunk <- 256L
  for (lo in seq.int(1L,ncol(mat),by=chunk)) {
    hi <- min(ncol(mat),lo+chunk-1L)
    writeBin(as.vector(as.matrix(mat[,lo:hi,drop=FALSE])),con,size=4,endian='little')
  }
  close(con);on.exit(NULL)
  stopifnot(file.rename(paste0(path,'.part'),path))
}
manifest_path <- file.path(dest,'manifest.json')
if (file.exists(manifest_path)) {
  prior <- fromJSON(manifest_path,simplifyVector=FALSE)
  stopifnot(identical(prior$input_sha256,input_sha),identical(prior$protocol_sha256,protocol_sha))
}
manifest <- list(status='running',group=group,job=Sys.getenv('SLURM_JOB_ID'),
  input_sha256=input_sha,protocol_sha256=protocol_sha,R=R.version.string,
  Seurat=as.character(packageVersion('Seurat')),SeuratObject=as.character(packageVersion('SeuratObject')),
  seed=42,norm='LogNormalize10000 over original RNA counts before geometry filtering',
  fixed_cells=TRUE,annotation_fields_in_fit=FALSE,samples=list())
write_manifest <- function() write_json(manifest,manifest_path,pretty=TRUE,auto_unbox=TRUE,na='null')
write_manifest()
payload <- readRDS(source_path)
counts <- payload$counts;qc <- payload$qc_metadata
stopifnot(identical(colnames(counts),rownames(qc)),!anyDuplicated(rownames(qc)))
stopifnot(all(Matrix::colSums(counts)==qc$nCount_RNA))
samples <- unique(qc$sample_id)
manifest$n_cells <- ncol(counts)
manifest$sample_order <- samples

# Shared full-RNA annotation input is normalized identically in pooled and single-cell-sample arms.
rna_path <- file.path(dest,'RNA_normalized.rds')
if (file.exists(rna_path)) {
  rna <- readRDS(rna_path)
} else {
  log_step('Creating full RNA annotation object')
  rna <- CreateSeuratObject(counts=counts,meta.data=qc,min.cells=0,min.features=0,project=group)
  rna <- NormalizeData(rna,normalization.method='LogNormalize',scale.factor=10000,verbose=FALSE)
  retained <- Matrix::rowSums(counts>0)>=3L
  rna <- subset(rna,features=rownames(counts)[retained])
  save_checkpoint(rna,rna_path)
}
rm(counts,payload);gc()
manifest$n_RNA_genes <- nrow(rna)
write.csv(rna@meta.data,file.path(dest,'cells.csv'))
writeLines(rownames(rna),file.path(dest,'RNA_genes.txt'))
objects <- list()
for (sample in samples) {
  sdir <- file.path(base,'samples',sample);dir.create(sdir,recursive=TRUE,showWarnings=FALSE)
  object_path <- file.path(sdir,'RNA_normalized_VST.rds')
  done_path <- file.path(sdir,'PREPARED')
  if (file.exists(done_path)) {
    stopifnot(identical(readLines(done_path),input_sha))
    obj <- readRDS(object_path)
  } else {
    log_step(paste('Preparing sample',sample))
    obj <- subset(rna,cells=rownames(rna@meta.data)[rna$sample_id==sample])
    detected <- Matrix::rowSums(obj[['RNA']]@counts>0)>=3L
    geometry_genes <- rownames(obj)[detected]
    # Fit VST on sample-eligible genes; annotation object retains all group RNA genes.
    geometry <- subset(obj,features=geometry_genes)
    geometry <- FindVariableFeatures(geometry,selection.method='vst',nfeatures=5000,verbose=FALSE)
    hvg <- VariableFeatures(geometry)
    feature_meta <- geometry[['RNA']]@meta.features
    feature_meta$gene <- rownames(feature_meta)
    feature_meta$hvg_rank_top5000 <- match(rownames(feature_meta),hvg)
    write.csv(feature_meta,file.path(sdir,'VST_feature_statistics.csv'),row.names=FALSE)
    writeLines(hvg,file.path(sdir,'HVG_rank_top5000.txt'))
    writeLines(geometry_genes,file.path(sdir,'geometry_genes.txt'))
    write.csv(obj@meta.data,file.path(sdir,'cells.csv'))
    log_step(paste('Scaling all geometry genes',sample,length(geometry_genes)))
    geometry <- ScaleData(geometry,features=geometry_genes,do.center=TRUE,do.scale=TRUE,
                          scale.max=10,verbose=FALSE)
    stopifnot(identical(rownames(geometry[['RNA']]@scale.data),geometry_genes))
    write_binary(geometry[['RNA']]@scale.data,file.path(sdir,'geometry_all.float32.bin'))
    write_json(list(n_cells=ncol(geometry),n_genes=length(geometry_genes),dtype='float32',
       order='C_cells_by_genes',endian='little',group=group,source_input_sha256=input_sha),
       file.path(sdir,'geometry_all.json'),auto_unbox=TRUE,pretty=TRUE)
    VariableFeatures(obj) <- hvg
    save_checkpoint(obj,object_path)
    writeLines(input_sha,done_path)
    rm(geometry);gc()
  }
  # The integration list retains only per-sample eligible genes, as in the reference workflow.
  geometry_genes <- readLines(file.path(sdir,'geometry_genes.txt'))
  integration_obj <- subset(obj,features=geometry_genes)
  VariableFeatures(integration_obj) <- head(readLines(file.path(sdir,'HVG_rank_top5000.txt')),2000L)
  objects[[sample]] <- integration_obj
  manifest$samples[[sample]] <- list(n_cells=ncol(obj),n_geometry_genes=length(geometry_genes),
                                    n_annotation_genes=nrow(obj))
  rm(obj,integration_obj);gc();write_manifest()
}
common_path <- file.path(dest,'integration_features.txt')
if (file.exists(common_path)) features <- readLines(common_path) else {
  log_step('Selecting 2000 shared integration features')
  features <- SelectIntegrationFeatures(object.list=objects,nfeatures=2000)
  writeLines(features,common_path)
}
stopifnot(length(features)==2000L,all(features %in% rownames(rna)))
manifest$integration_feature_count <- length(features)
write_binary(rna[['RNA']]@data[features,,drop=FALSE],file.path(dest,'DL_RNA2000.float32.bin'))
write_json(list(n_cells=ncol(rna),n_genes=length(features),dtype='float32',order='C_cells_by_genes',
  endian='little',gene_file='integration_features.txt',cell_file='cells.csv'),
  file.path(dest,'DL_RNA2000.json'),auto_unbox=TRUE,pretty=TRUE)

none_path <- file.path(dest,'NONE_PCA30.rds')
if (!file.exists(none_path)) {
  log_step('Computing matched uncorrected RNA PCA30')
  VariableFeatures(rna) <- features
  rna <- ScaleData(rna,features=features,do.center=TRUE,do.scale=TRUE,scale.max=10,verbose=FALSE)
  rna <- RunPCA(rna,features=features,npcs=30,seed.use=42,verbose=FALSE)
  save_checkpoint(rna,none_path)
  write.csv(Embeddings(rna,'pca'),file.path(dest,'NONE_PCA30.csv'))
} else if (!file.exists(file.path(dest,'NONE_PCA30.csv'))) {
  recovered <- readRDS(none_path)
  write.csv(Embeddings(recovered,'pca'),file.path(dest,'NONE_PCA30.csv'))
  rm(recovered);gc()
}
rm(rna);gc()
anchor_path <- file.path(dest,'CCA_anchors.rds')
if (file.exists(anchor_path)) anchors <- readRDS(anchor_path) else {
  log_step('Finding CCA anchors across the four samples')
  anchors <- FindIntegrationAnchors(object.list=objects,anchor.features=features,
      normalization.method='LogNormalize',reduction='cca',dims=1:30,
      k.anchor=5,k.filter=200,k.score=30,max.features=200,nn.method='annoy',n.trees=50,
      scale=TRUE,verbose=TRUE)
  save_checkpoint(anchors,anchor_path)
}
rm(objects);gc()
cca_path <- file.path(dest,'CCA_PCA30.rds')
if (!file.exists(cca_path)) {
  log_step('Integrating CCA expression and computing PCA30')
  integrated <- IntegrateData(anchorset=anchors,new.assay.name='integrated',
      normalization.method='LogNormalize',dims=1:30,k.weight=100,sd.weight=1,
      preserve.order=FALSE,verbose=TRUE)
  DefaultAssay(integrated) <- 'integrated'
  integrated <- ScaleData(integrated,features=features,do.center=TRUE,do.scale=TRUE,
                           scale.max=10,verbose=FALSE)
  integrated <- RunPCA(integrated,features=features,npcs=30,seed.use=42,verbose=FALSE)
  stopifnot(setequal(colnames(integrated),rownames(qc)))
  save_checkpoint(integrated,cca_path)
  write.csv(Embeddings(integrated,'pca'),file.path(dest,'CCA_PCA30.csv'))
  # Explicit cell ordering matches the RNA input; corrected input is a separate legacy arm.
  write_binary(integrated[['integrated']]@data[features,rownames(qc),drop=FALSE],
               file.path(dest,'DL_CCA2000.float32.bin'))
} else if (!file.exists(file.path(dest,'DL_CCA2000.float32.bin')) ||
           !file.exists(file.path(dest,'CCA_PCA30.csv'))) {
  integrated <- readRDS(cca_path)
  write.csv(Embeddings(integrated,'pca'),file.path(dest,'CCA_PCA30.csv'))
  write_binary(integrated[['integrated']]@data[features,rownames(qc),drop=FALSE],
               file.path(dest,'DL_CCA2000.float32.bin'))
}
manifest$status <- 'completed'
manifest$elapsed_seconds <- as.numeric(difftime(Sys.time(),start,units='secs'))
manifest$completed_at <- format(Sys.time(),tz='UTC')
write_manifest()
writeLines(capture.output(sessionInfo()),file.path(dest,'sessionInfo.txt'))
writeLines(digest(file=manifest_path,algo='sha256'),file.path(dest,'PREPARATION_COMPLETE'))
log_step('Fixed-QC preprocessing and CCA preparation complete')
