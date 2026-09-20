stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
suppressPackageStartupMessages(library(Seurat));suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest));suppressPackageStartupMessages(library(future))
plan(sequential)
root<-'/fs/scratch/PCON0080/yimin/dgscrna';base<-file.path(root,'results/hvg_ptc_20260916_v1/paper_claim_validation_20260917')
a<-commandArgs(trailingOnly=TRUE);cfg<-fromJSON(a[[1]])
source<-file.path(base,'GBM',cfg$sample,cfg$budget)
dest<-cfg$dest
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
if(file.exists(file.path(dest,'PREPARED')))quit(status=0)
m<-fromJSON(file.path(source,'prepare_manifest.json'))
stopifnot(readLines(file.path(source,'PREPARED'))==digest(file=file.path(source,'prepare_manifest.json'),algo='sha256'))
obj<-readRDS(file.path(source,'expression_PCA30.rds'))
if(cfg$embedding_seed!=42L) {
  # Feature identities and normalized DL stay fixed. This perturbation reruns
  # PCA and (when used) UMAP; the separate MLP seed remains42.
  features<-readLines(file.path(source,'selected_features.txt'))
  obj<-RunPCA(obj,features=features,npcs=30L,seed.use=as.integer(cfg$embedding_seed),verbose=FALSE)
}
if(startsWith(cfg$space,'UMAP')) {
  dims<-as.integer(sub('UMAP','',cfg$space))
  if(dims!=2L || cfg$embedding_seed!=42L || cfg$umap_neighbors!=30L)obj<-RunUMAP(obj,reduction='pca',dims=1:30,
    n.components=dims,n.neighbors=as.integer(cfg$umap_neighbors),umap.method='uwot',metric='cosine',min.dist=.3,
    seed.use=as.integer(cfg$embedding_seed),verbose=FALSE)
  Z<-Embeddings(obj,'umap')
} else if(cfg$space=='PCA30')Z<-Embeddings(obj,'pca') else if(cfg$space=='RNA_noDR') {
  Z<-t(as.matrix(obj[['RNA']]@scale.data))
} else stop('Unsupported frozen control space')
stopifnot(all(is.finite(Z)),identical(rownames(Z),colnames(obj)))
writeLines(if(cfg$space=='RNA_noDR')rownames(obj[['RNA']]@scale.data) else rownames(Loadings(obj,'pca')),
  file.path(dest,'control_feature_genes.txt'))
colnames(Z)<-paste0('CTRL_',seq_len(ncol(Z)))
obj[['control']]<-CreateDimReducObject(embeddings=Z,key='CTRL_',assay='RNA')
write.csv(Z,file.path(dest,'control_embedding.csv'))
saveRDS(obj,file.path(dest,'expression_PCA30.rds'),compress=FALSE)
for(file in c('cells.csv','selected_features.txt','geometry_features.txt','scoring_features.txt','DL_features.txt')) {
  target<-file.path(dest,file)
  if(!file.exists(target))stopifnot(file.link(file.path(source,file),target))
  stopifnot(digest(file=target,algo='sha256')==digest(file=file.path(source,file),algo='sha256'))
}
m$reference_prepare_manifest_sha256<-digest(file=file.path(source,'prepare_manifest.json'),algo='sha256')
m$expression_sha256<-digest(file=file.path(dest,'expression_PCA30.rds'),algo='sha256')
m$representation_control<-cfg;m$geometry_dimensions<-ncol(Z)
m$DL_binary<-file.path(source,'DL.float32.bin');m$seed<-cfg$embedding_seed
m$job<-Sys.getenv('SLURM_JOB_ID');m$condition<-cfg$name
m$scientific_scope<-'Fixed feature/scoring/DL inputs within a budget; explicitly varied embedding dimension, clustering parameter, or PCA/UMAP seed.'
write_json(m,file.path(dest,'prepare_manifest.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(digest(file=file.path(dest,'prepare_manifest.json'),algo='sha256'),file.path(dest,'PREPARED'))
