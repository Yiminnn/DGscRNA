#!/usr/bin/env Rscript
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest))
suppressPackageStartupMessages(library(future))
plan(sequential)
options(Seurat.object.assay.version='v3',future.globals.maxSize=64*1024^3)
root<-'/fs/scratch/PCON0080/yimin/dgscrna'
base<-file.path(root,'results/hvg_ptc_20260916_v1/paper_claim_validation_20260917')
args<-commandArgs(trailingOnly=TRUE)
sample<-args[[1]];budget<-args[[2]]
stopifnot(budget %in% c('hvg500','hvg1000','hvg2000','hvg3000','hvg5000','all'))
seed<-if(length(args)>2L)as.integer(args[[3]]) else 42L
condition<-if(seed==42L)budget else paste0(budget,'_seed',seed)
inp<-file.path(base,'inputs',sample)
dest<-file.path(base,'GBM',sample,condition);dir.create(dest,recursive=TRUE,showWarnings=FALSE)
im<-fromJSON(file.path(inp,'input_manifest.json'))
stopifnot(readLines(file.path(inp,'INPUT_COMPLETE'))==digest(file=file.path(inp,'input_manifest.json'),algo='sha256'))
if(file.exists(file.path(dest,'PREPARED'))) {
  stopifnot(readLines(file.path(dest,'PREPARED'))==digest(file=file.path(dest,'prepare_manifest.json'),algo='sha256'))
  quit(status=0)
}
started<-Sys.time();set.seed(seed)
log<-function(...) {cat(format(Sys.time()),..., '\n');flush.console()}
publish<-function(path,writer) {tmp<-paste0(path,'.part.',Sys.getenv('SLURM_JOB_ID'));writer(tmp);stopifnot(file.rename(tmp,path))}
save_atomic<-function(x,path) publish(path,function(p)saveRDS(x,p,compress=FALSE))
cells<-read.csv(file.path(inp,'cells_fit.csv'),stringsAsFactors=FALSE,check.names=FALSE)
genes<-read.csv(file.path(inp,'genes.csv'),stringsAsFactors=FALSE)$gene
stopifnot(!anyDuplicated(cells$cell_id),!anyDuplicated(genes),length(unique(cells$batch))==1L)
read_bin<-function(name,what,n,size) {con<-file(file.path(inp,name),'rb');on.exit(close(con));readBin(con,what,n=n,size=size,endian='little')}
for(n in names(im$fitting_files))stopifnot(digest(file=file.path(inp,n),algo='sha256')==im$fitting_files[[n]])
counts<-new('dgCMatrix',x=read_bin('x.bin','numeric',im$nnz,8),i=read_bin('i.bin','integer',im$nnz,4),
  p=read_bin('p.bin','integer',im$n_cells+1L,4),Dim=as.integer(c(im$n_genes,im$n_cells)),Dimnames=list(genes,cells$cell_id))
rownames(cells)<-cells$cell_id
obj<-CreateSeuratObject(counts,meta.data=cells,min.cells=0,min.features=0,project=sample)
stopifnot(nrow(obj)==length(genes),!anyDuplicated(rownames(obj)))
publish(file.path(dest,'Seurat_gene_names.csv'),function(p)write.csv(data.frame(source=genes,Seurat=rownames(obj)),p,row.names=FALSE))
obj<-NormalizeData(obj,normalization.method='LogNormalize',scale.factor=10000,verbose=FALSE)
rm(counts);gc()
requested<-if(budget=='all')nrow(obj) else as.integer(sub('hvg','',budget))
if(budget=='all') {
  features<-rownames(obj)
} else {
  obj<-FindVariableFeatures(obj,selection.method='vst',nfeatures=requested,loess.span=0.3,clip.max='auto',verbose=FALSE)
  features<-VariableFeatures(obj)
  publish(file.path(dest,'native_vst_statistics.csv'),function(p)write.csv(obj[['RNA']]@meta.features,p))
}
stopifnot(length(features)>=31L,ncol(obj)>31L)
VariableFeatures(obj)<-features
log('SCALE_PCA',sample,budget,length(features),ncol(obj))
obj<-ScaleData(obj,features=features,do.center=TRUE,do.scale=TRUE,scale.max=10,verbose=FALSE)
obj<-RunPCA(obj,features=features,npcs=30,seed.use=seed,verbose=FALSE)
stopifnot(ncol(Embeddings(obj,'pca'))==30L)
obj<-RunUMAP(obj,reduction='pca',dims=1:30,n.components=2,n.neighbors=30L,
  umap.method='uwot',metric='cosine',min.dist=0.3,seed.use=seed,verbose=FALSE)
actual_geometry<-rownames(Loadings(obj,'pca'))
writeLines(features,file.path(dest,'selected_features.txt'))
writeLines(actual_geometry,file.path(dest,'geometry_features.txt'))
# Native R single-batch branch scores every RNA gene and trains DL on selected
# normalized genes. Do not incorrectly claim HVG also limits RNA DEG scoring.
features<-rownames(obj[['RNA']]@data)[rownames(obj[['RNA']]@data) %in% features]
writeLines(features,file.path(dest,'DL_features.txt'))
writeLines(rownames(obj[['RNA']]@data),file.path(dest,'scoring_features.txt'))
expression<-obj[['RNA']]@data[features,,drop=FALSE]
publish(file.path(dest,'DL.float32.bin'),function(p) {
  con<-file(p,'wb');on.exit(close(con))
  for(lo in seq.int(1L,ncol(expression),by=512L)) {
    hi<-min(ncol(expression),lo+511L)
    writeBin(as.numeric(as.matrix(expression[,lo:hi,drop=FALSE])),con,size=4,endian='little')
  }
})
rm(expression);gc()
publish(file.path(dest,'cells.csv'),function(p)write.csv(data.frame(cell_id=colnames(obj),batch=obj$batch),p,row.names=FALSE))
publish(file.path(dest,'PCA30.csv'),function(p)write.csv(Embeddings(obj,'pca'),p))
publish(file.path(dest,'UMAP2.csv'),function(p)write.csv(Embeddings(obj,'umap'),p))
save_atomic(obj,file.path(dest,'expression_PCA30.rds'))
script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[[1]])
m<-list(status='completed',unit=sample,dataset='GSE274546',sample=sample,condition=condition,budget=budget,
  n_cells=ncol(obj),seed=seed,requested_features=requested,
  features=list(anchor=0L,geometry=length(actual_geometry),selected=length(features),scoring=nrow(obj[['RNA']]),DL=length(features)),
  correction='no_correction_single_batch',assay='RNA',n_batches=1L,
  input_manifest_sha256=digest(file=file.path(inp,'input_manifest.json'),algo='sha256'),input_semantics=im$input_semantics,
  DL_binary=file.path(dest,'DL.float32.bin'),DL_binary_sha256=digest(file=file.path(dest,'DL.float32.bin'),algo='sha256'),
  expression_sha256=digest(file=file.path(dest,'expression_PCA30.rds'),algo='sha256'),
  stage_gene_sha256=lapply(c('selected_features.txt','geometry_features.txt','scoring_features.txt','DL_features.txt'),function(n)list(file=n,sha256=digest(file=file.path(dest,n),algo='sha256'))),
  scientific_scope='Native R RNA single-sample HVG changes geometry and DL; scoring remains all RNA genes. No CCA across patients.',
  QC='Existing author cells and >=3-cell gene filter; no additional selection',
  reference_labels_used_for_fitting=FALSE,execution_source=script,source_sha256=digest(file=script,algo='sha256'),
  job=Sys.getenv('SLURM_JOB_ID'),elapsed_seconds=as.numeric(difftime(Sys.time(),started,units='secs')))
publish(file.path(dest,'prepare_manifest.json'),function(p)write_json(m,p,pretty=TRUE,auto_unbox=TRUE))
writeLines(capture.output(sessionInfo()),file.path(dest,'sessionInfo.txt'))
publish(file.path(dest,'PREPARED'),function(p)writeLines(digest(file=file.path(dest,'prepare_manifest.json'),algo='sha256'),p))
log('PREPARED',sample,budget)
