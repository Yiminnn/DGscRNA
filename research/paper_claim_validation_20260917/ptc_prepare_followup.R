#!/usr/bin/env Rscript
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
suppressPackageStartupMessages(library(jsonlite));suppressPackageStartupMessages(library(digest))
root<-'/fs/scratch/PCON0080/yimin/dgscrna'
base<-file.path(root,'results/hvg_ptc_20260916_v1/paper_claim_validation_20260917')
old<-file.path(root,'results/hvg_ptc_20260916_v1/r_reference_campaign_20260917')
gate<-fromJSON(file.path(base,'GBM_full_summary/GBM_FULL_DELIVERED.json'))
stopifnot(isTRUE(gate$PTC_compute_may_start),gate$receipt_sha256==digest(file=file.path(base,'GBM_full_summary/DELIVERY_RECEIPT.json'),algo='sha256'))
suppressPackageStartupMessages(library(Seurat));suppressPackageStartupMessages(library(future))
plan(sequential);options(Seurat.object.assay.version='v3')
cfg<-fromJSON(commandArgs(trailingOnly=TRUE)[[1]])
dest<-cfg$dest;dir.create(dest,recursive=TRUE,showWarnings=FALSE)
if(file.exists(file.path(dest,'PREPARED'))) {
  stopifnot(readLines(file.path(dest,'PREPARED'))==digest(file=file.path(dest,'prepare_manifest.json'),algo='sha256'))
  quit(status=0)
}
stopifnot(cfg$group %in% c('NMT','TTU'),cfg$kind %in% c('representation','retention','default_parity'))
verified_manifest<-function(p) {
  stopifnot(readLines(file.path(p,'PREPARED'))==digest(file=file.path(p,'prepare_manifest.json'),algo='sha256'))
  fromJSON(file.path(p,'prepare_manifest.json'))
}
link_exact<-function(a,b) {
  if(!file.exists(b))stopifnot(file.link(a,b))
  stopifnot(digest(file=a,algo='sha256')==digest(file=b,algo='sha256'))
}
if(cfg$kind=='retention') {
  source<-file.path(old,'PTC_ablation',paste0('PTC_',cfg$group,'_CCAall'))
  pm<-verified_manifest(source)
  obj<-readRDS(file.path(source,'expression_PCA30.rds'))
  all_features<-readLines(file.path(source,'DL_features.txt'))
  fixed_dir<-file.path(old,'PTC_geometry_fixed_inputs',cfg$group)
  fixed<-readLines(file.path(fixed_dir,'DL_features.txt'))
  marker_path<-file.path(old,'markers/PTC_original17.json')
  libs<-fromJSON(marker_path,simplifyVector=FALSE)
  marker_genes<-unique(unlist(libs,use.names=FALSE))
  stopifnot(DefaultAssay(obj)=='integrated',identical(all_features,rownames(obj[['integrated']])) )
  # Keep fixed2000's original order, then append eligible extra markers in CCAall order.
  # This also preserves common genes' relative tie/summation order in the old scorer.
  scoring<-c(fixed,all_features[all_features %in% setdiff(marker_genes,fixed)])
  stopifnot(all(fixed %in% scoring),length(fixed)==2000L)
  assay<-CreateAssayObject(data=obj[['integrated']]@data[scoring,,drop=FALSE])
  result<-CreateSeuratObject(counts=assay,assay='CCAall_marker_retained',meta.data=obj@meta.data)
  rm(obj,assay);gc()
  geometry<-file.path(old,'PTC_ablation',paste0('PTC_',cfg$group,'_GEOMETRY',cfg$budget,'_FIXED_CCAall_DL2000'))
  gm<-verified_manifest(geometry)
  stopifnot(gm$fixed_expression_source_sha256==digest(file=file.path(source,'prepare_manifest.json'),algo='sha256'))
  pca<-as.matrix(read.csv(file.path(geometry,'PCA30.csv'),row.names=1,check.names=FALSE))
  umap<-as.matrix(read.csv(file.path(geometry,'UMAP2.csv'),row.names=1,check.names=FALSE))
  stopifnot(identical(rownames(pca),colnames(result)),identical(rownames(umap),colnames(result)))
  result[['pca']]<-CreateDimReducObject(embeddings=pca,key='PC_',assay=DefaultAssay(result))
  result[['umap']]<-CreateDimReducObject(embeddings=umap,key='UMAP_',assay=DefaultAssay(result))
  # Exact existing partitions make this a scoring-universe, not clustering, perturbation.
  for(route in c('PCA30_SNN','UMAP2_HDBSCAN_R')) {
    rd<-file.path(dest,route);dir.create(rd,showWarnings=FALSE)
    for(f in c('clusters.csv','density_diagnostics.csv'))
      if(file.exists(file.path(geometry,route,f)))link_exact(file.path(geometry,route,f),file.path(rd,f))
  }
  for(f in c('PCA30.csv','UMAP2.csv','cells.csv','geometry_features.txt'))link_exact(file.path(geometry,f),file.path(dest,f))
  link_exact(file.path(fixed_dir,'DL_features.txt'),file.path(dest,'DL_features.txt'))
  pm$DL_binary<-file.path(fixed_dir,'DL.float32.bin');pm$DL_binary_sha256<-digest(file=pm$DL_binary,algo='sha256')
  pm$features<-list(anchor=length(all_features),geometry=gm$features$geometry,scoring=length(scoring),DL=2000L)
  pm$correction<-'fixed_CCAall_uniform_marker_retention_DL2000'
  pm$marker_source_sha256<-digest(file=marker_path,algo='sha256')
  pm$geometry_source<-geometry;pm$geometry_source_manifest_sha256<-digest(file=file.path(geometry,'prepare_manifest.json'),algo='sha256')
  writeLines(scoring,file.path(dest,'scoring_features.txt'))
  writeLines(sort(setdiff(marker_genes,all_features)),file.path(dest,'markers_outside_eligible_CCA.txt'))
  pm$uniform_retention_rule<-'(fixed2000 union all17 marker genes) intersect eligible group CCAall universe; fixed2000 order then extras in CCAall order'
  pm$scientific_scope<-'Conditional on frozen all-gene CCA expression; only score-gene universe changed from original fixed2000 geometry control'
} else {
  source<-file.path(old,'PTC_ablation',paste0('PTC_',cfg$group,'_CCA',cfg$budget))
  pm<-verified_manifest(source);result<-readRDS(file.path(source,'expression_PCA30.rds'))
  features<-readLines(file.path(source,'DL_features.txt'))
  pca_features<-rownames(Loadings(result,'pca'))
  stopifnot(length(pca_features)>=30L,all(pca_features %in% features))
  if(cfg$seed==42L) {
    umap<-as.matrix(read.csv(file.path(source,'UMAP2.csv'),row.names=1,check.names=FALSE))
    stopifnot(identical(rownames(umap),colnames(result)))
    result[['umap']]<-CreateDimReducObject(embeddings=umap,key='UMAP_',assay=DefaultAssay(result))
  } else {
    # Preserve the actually used original PCA feature order, including its zero-variance exclusions.
    result<-RunPCA(result,features=pca_features,npcs=30L,seed.use=as.integer(cfg$seed),verbose=FALSE)
    result<-RunUMAP(result,reduction='pca',dims=1:30,n.components=2,n.neighbors=30L,
      umap.method='uwot',metric='cosine',min.dist=.3,seed.use=as.integer(cfg$seed),verbose=FALSE)
  }
  for(f in c('cells.csv','DL_features.txt'))link_exact(file.path(source,f),file.path(dest,f))
  writeLines(features,file.path(dest,'geometry_features.txt'));writeLines(rownames(result[[DefaultAssay(result)]]),file.path(dest,'scoring_features.txt'))
  writeLines(pca_features,file.path(dest,'PCA_actual_feature_order.txt'))
  write.csv(Embeddings(result,'pca'),file.path(dest,'PCA30.csv'))
  write.csv(Embeddings(result,'umap'),file.path(dest,'UMAP2.csv'))
  pm$scientific_scope<-'PCA/UMAP seed sensitivity conditional on unchanged CCA fit, features and expression; SNN seed0 and MLP seed42 fixed'
}
cells<-read.csv(file.path(dest,'cells.csv'),stringsAsFactors=FALSE)
stopifnot(identical(cells$cell_id,colnames(result)),all(is.finite(Embeddings(result,'pca'))),all(is.finite(Embeddings(result,'umap'))))
saveRDS(result,file.path(dest,'expression_PCA30.rds.part'),compress=FALSE)
stopifnot(file.rename(file.path(dest,'expression_PCA30.rds.part'),file.path(dest,'expression_PCA30.rds')))
pm$unit<-cfg$name;pm$condition<-cfg;pm$assay<-DefaultAssay(result);pm$seed<-cfg$seed
pm$source_preparation<-source;pm$source_prepare_manifest_sha256<-digest(file=file.path(source,'prepare_manifest.json'),algo='sha256')
pm$expression_sha256<-digest(file=file.path(dest,'expression_PCA30.rds'),algo='sha256')
pm$job<-Sys.getenv('SLURM_JOB_ID');pm$reference_labels_used_for_fitting<-FALSE
pm$execution_source<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[[1]])
pm$source_sha256<-digest(file=pm$execution_source,algo='sha256')
write_json(pm,file.path(dest,'prepare_manifest.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),file.path(dest,'sessionInfo.txt'))
writeLines(digest(file=file.path(dest,'prepare_manifest.json'),algo='sha256'),file.path(dest,'PREPARED'))
