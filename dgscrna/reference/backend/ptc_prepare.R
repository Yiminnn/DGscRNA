#!/usr/bin/env Rscript
# Advanced prepared-input CCA backend. This does not perform historical cell QC.
if (Sys.getenv('DGSCRNA_REQUIRE_SLURM','0')=='1' && !nzchar(Sys.getenv('SLURM_JOB_ID')))
  stop('This execution requires a SLURM allocation')
execution_id <- function() {
  id<-Sys.getenv('SLURM_JOB_ID')
  if(nzchar(id)) id else paste0('local-',Sys.getpid())
}
script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[[1]])
reference_lib<-Sys.getenv('DGSCRNA_REFERENCE_R_LIB','')
reference_libraries<-c(if(nzchar(reference_lib))reference_lib else character(),.Library)
.libPaths(reference_libraries,include.site=FALSE)
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest))
suppressPackageStartupMessages(library(future))
plan(sequential);options(Seurat.object.assay.version='v3',future.globals.maxSize=250*1024^3)
base<-Sys.getenv('DGSCRNA_REFERENCE_OUT',Sys.getenv('DGSCRNA_EXAMPLE_OUT'));stopifnot(nzchar(base))
args<-commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L)
config_path<-normalizePath(args[[1]],mustWork=TRUE)
spec<-fromJSON(config_path,simplifyVector=FALSE)
stopifnot(spec$group %in% c('NMT','TTU'),spec$correction=='CCA')
sample_groups<-list(NMT=c('MT-1','MT-2','N-1','N-2'),TTU=c('TU-1','TU-2','T-1','T-2'))
samples<-sample_groups[[spec$group]]
stopifnot(identical(names(spec$samples),samples))
unit<-if(is.null(spec$unit))spec$group else spec$unit
stopifnot(grepl('^[A-Za-z0-9][A-Za-z0-9_.-]*$',unit))
budget<-if(spec$hvg=='all')'all' else paste0('hvg',spec$hvg)
stopifnot(budget %in% c('hvg500','hvg1000','hvg2000','hvg3000','hvg5000','all'))
seed<-if(is.null(spec$seed))42L else as.integer(spec$seed)
stopifnot(length(seed)==1L,!is.na(seed),seed>=0L)
condition<-if(seed==42L)budget else paste0(budget,'_seed',seed)
dest<-file.path(base,'PTC',unit,condition);dir.create(dest,recursive=TRUE,showWarnings=FALSE)
resolve_input<-function(path) {
  if(!grepl('^/',path))path<-file.path(dirname(config_path),path)
  normalizePath(path,mustWork=TRUE)
}
input_paths<-list(canonical_cells=resolve_input(spec$canonical_cells))
for(sample in samples) {
  input_paths[[paste0(sample,'_object')]]<-resolve_input(spec$samples[[sample]]$normalized_rds)
  input_paths[[paste0(sample,'_genes')]]<-resolve_input(spec$samples[[sample]]$geometry_genes)
}
input_hashes<-lapply(input_paths,function(p)digest(file=p,algo='sha256'))
config_sha256<-digest(file=config_path,algo='sha256')
if(file.exists(file.path(dest,'PREPARED'))) {
  stopifnot(readLines(file.path(dest,'PREPARED'))==digest(file=file.path(dest,'prepare_manifest.json'),algo='sha256'))
  prior<-fromJSON(file.path(dest,'prepare_manifest.json'),simplifyVector=FALSE)
  stopifnot(identical(prior$prepared_input_sha256,input_hashes),prior$config_sha256==config_sha256)
  quit(status=0)
}
started<-Sys.time();set.seed(seed)
publish<-function(path,writer) {
  temporary<-paste0(path,'.part.',execution_id(),'.',Sys.getpid())
  writer(temporary);stopifnot(file.rename(temporary,path))
}
save_atomic<-function(x,path)publish(path,function(p)saveRDS(x,p,compress=FALSE))
objects<-list()
for(sample in samples) {
  o<-readRDS(input_paths[[paste0(sample,'_object')]])
  eligible<-readLines(input_paths[[paste0(sample,'_genes')]])
  stopifnot(inherits(o,'Seurat'),!anyDuplicated(eligible),all(eligible %in% rownames(o[['RNA']])) )
  DefaultAssay(o)<-'RNA'
  stopifnot(ncol(o)>100L,ncol(o[['RNA']]@data)==ncol(o))
  o<-subset(o,features=eligible)
  # Remove archived annotations from the fitting object, as in the source.
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
writeLines(features,file.path(dest,'selected_features.txt'))
# Anchor caches are only accepted when all prepared inputs and this configuration match.
anchorpath<-file.path(dest,'anchors.rds')
anchor_manifest_path<-file.path(dest,'anchors_manifest.json')
if(file.exists(anchorpath) && file.exists(anchor_manifest_path)) {
  am<-fromJSON(anchor_manifest_path,simplifyVector=FALSE)
  stopifnot(identical(am$prepared_input_sha256,input_hashes),am$config_sha256==config_sha256,
            am$anchors_sha256==digest(file=anchorpath,algo='sha256'))
  anchors<-readRDS(anchorpath)
} else {
  anchors<-FindIntegrationAnchors(objects,anchor.features=features,normalization.method='LogNormalize',reduction='cca',
     dims=1:30,k.anchor=5,k.filter=200,k.score=30,max.features=200,nn.method='annoy',n.trees=50,scale=TRUE,verbose=TRUE)
  save_atomic(anchors,anchorpath)
  publish(anchor_manifest_path,function(p)write_json(list(prepared_input_sha256=input_hashes,
    config_sha256=config_sha256,anchors_sha256=digest(file=anchorpath,algo='sha256')),p,auto_unbox=TRUE,pretty=TRUE))
}
rm(objects);gc()
obj<-IntegrateData(anchors,new.assay.name='integrated',normalization.method='LogNormalize',dims=1:30,
      k.weight=100,sd.weight=1,verbose=TRUE)
rm(anchors);gc();DefaultAssay(obj)<-'integrated'
obj<-ScaleData(obj,features=features,do.center=TRUE,do.scale=TRUE,scale.max=10,verbose=FALSE)
obj<-RunPCA(obj,features=features,npcs=30,seed.use=seed,verbose=FALSE)
correction<-paste0('fresh_Seurat_CCA_',spec$group,'_',spec$hvg)
canonical<-read.csv(input_paths$canonical_cells,stringsAsFactors=FALSE,check.names=FALSE,
  colClasses='character',na.strings=NULL)
stopifnot('cell_id' %in% names(canonical))
canonical_cells<-canonical$cell_id
stopifnot(!anyDuplicated(canonical_cells),setequal(canonical_cells,colnames(obj)))
order_changed<-!identical(canonical_cells,colnames(obj))
if(order_changed) {
  obj<-obj[,canonical_cells]
  stopifnot(identical(colnames(obj),canonical_cells))
  obj<-RunPCA(obj,features=features,npcs=30,seed.use=seed,verbose=FALSE)
}
assay<-DefaultAssay(obj)
features<-rownames(obj[[assay]]@data)[rownames(obj[[assay]]@data) %in% features]
writeLines(features,file.path(dest,'DL_features.txt'))
writeLines(rownames(Loadings(obj,'pca')),file.path(dest,'geometry_features.txt'))
writeLines(rownames(obj[[assay]]@data),file.path(dest,'scoring_features.txt'))
writeLines(correction,file.path(dest,'correction.txt'))
prep_path<-file.path(dest,'expression_PCA30.rds')
save_atomic(obj,prep_path)
publish(file.path(dest,'DL.float32.bin'),function(p) {
  con<-file(p,'wb');on.exit(close(con))
  for(lo in seq.int(1L,ncol(obj),by=256L)) {
    hi<-min(ncol(obj),lo+255L)
    writeBin(as.numeric(as.matrix(obj[[assay]]@data[features,lo:hi,drop=FALSE])),con,size=4,endian='little')
  }
})
publish(file.path(dest,'cells.csv'),function(p)write.csv(data.frame(cell_id=colnames(obj),batch=obj$sample_id),p,row.names=FALSE))
publish(file.path(dest,'PCA30.csv'),function(p)write.csv(Embeddings(obj,'pca'),p))
m<-list(status='completed',unit=unit,dataset='PTC',group=spec$group,condition=condition,budget=budget,
 n_cells=ncol(obj),correction=correction,n_batches=length(samples),
 features=list(anchor=length(features),geometry=length(rownames(Loadings(obj,'pca'))),selected=length(features),scoring=nrow(obj[[assay]]),DL=length(features)),
 assay=assay,seed=seed,DL_binary=file.path(dest,'DL.float32.bin'),DL_binary_sha256=digest(file=file.path(dest,'DL.float32.bin'),algo='sha256'),
 expression_sha256=digest(file=prep_path,algo='sha256'),
 prepared_input_sha256=input_hashes,config_sha256=config_sha256,prepared_input_paths=input_paths,
 reference_labels_used_for_fitting=FALSE,prepared_input_cells_not_refiltered=TRUE,historical_cell_QC_verified_by_this_stage=FALSE,
 canonical_cell_order_locked=TRUE,cell_order_changed_from_native_integration=order_changed,
 input_semantics='User-supplied prepared RNA Seurat objects normalized on full counts before per-sample geometry filtering; historical QC is external',
 no_HVG_definition=if(spec$hvg=='all')'all genes detected in >=3 cells in every included sample; no variance ranking' else NULL,
 QC='Prepared-input only; does not perform or independently verify historical cell QC',
 job=execution_id(),execution_source=script,source_sha256=digest(file=script,algo='sha256'),
 runtime_R_library_paths=.libPaths(),
 runtime_R_package_paths=setNames(lapply(c('Seurat','Matrix','jsonlite','digest','future'),function(p)normalizePath(find.package(p))),
                                 c('Seurat','Matrix','jsonlite','digest','future')),
 elapsed_seconds=as.numeric(difftime(Sys.time(),started,units='secs')))
publish(file.path(dest,'prepare_manifest.json'),function(p)write_json(m,p,auto_unbox=TRUE,pretty=TRUE))
writeLines(capture.output(sessionInfo()),file.path(dest,'sessionInfo.txt'))
publish(file.path(dest,'PREPARED'),function(p)writeLines(digest(file=file.path(dest,'prepare_manifest.json'),algo='sha256'),p))
cat('PREPARED',unit,ncol(obj),length(features),'\n')
