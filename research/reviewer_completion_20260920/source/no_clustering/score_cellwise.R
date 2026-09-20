#!/usr/bin/env Rscript
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest))
args<-commandArgs(trailingOnly=TRUE);cfg<-fromJSON(args[[1]],simplifyVector=FALSE)
prep<-cfg$prep;dest<-cfg$route_dir;dir.create(dest,recursive=TRUE,showWarnings=FALSE)
hash<-function(p)digest(file=p,algo='sha256')
write_atomic<-function(path,writer) {tmp<-paste0(path,'.part.',Sys.getpid());writer(tmp);stopifnot(file.rename(tmp,path))}
pm<-fromJSON(file.path(prep,'prepare_manifest.json'))
stopifnot(readLines(file.path(prep,'PREPARED'))==hash(file.path(prep,'prepare_manifest.json')),
          hash(file.path(prep,'prepare_manifest.json'))==cfg$input$prepare_manifest_sha256,
          hash(file.path(prep,'expression_PCA30.rds'))==cfg$input$expression_sha256,
          hash(file.path(prep,'cells.csv'))==cfg$input$cells_sha256,
          hash(cfg$marker_file)==cfg$marker_sha256)
if(file.exists(file.path(dest,'SCORE_COMPLETE'))) {
  stopifnot(readLines(file.path(dest,'SCORE_COMPLETE'))==hash(file.path(dest,'score_manifest.json')))
  old<-fromJSON(file.path(dest,'score_manifest.json'))
  stopifnot(old$input_signature==cfg$input_signature,old$initial_sha256==hash(file.path(dest,'initial_calls.csv.gz')))
  quit(status=0)
}
started<-Sys.time();obj<-readRDS(file.path(prep,'expression_PCA30.rds'))
cells<-read.csv(file.path(prep,'cells.csv'),stringsAsFactors=FALSE)
stopifnot(identical(colnames(obj),cells$cell_id),DefaultAssay(obj)=='RNA')
x<-obj[['RNA']]@data
stopifnot(inherits(x,'sparseMatrix'),all(is.finite(x@x)),all(x@x>=0),identical(colnames(x),cells$cell_id))
libraries<-fromJSON(cfg$marker_file,simplifyVector=FALSE)
panels<-lapply(libraries[['CM2_glioma_other']],unlist,use.names=FALSE)
stopifnot(length(panels)>1L,all(lengths(panels)>0L))
S<-matrix(0,length(panels),ncol(x),dimnames=list(names(panels),colnames(x)))
retention<-list()
for(p in names(panels)) {
  genes<-intersect(panels[[p]],rownames(x));denominator<-length(panels[[p]])
  if(length(genes))S[p,]<-Matrix::colSums(x[genes,,drop=FALSE])/denominator
  if(denominator<=1L)S[p,]<-S[p,]*0.8
  retention[[p]]<-data.frame(panel=p,denominator=denominator,retained=length(genes),
                            retained_genes=paste(genes,collapse=';'))
}
stopifnot(all(is.finite(S)),all(S>=0))
mx<-apply(S,2,max);ties<-colSums(S==rep(mx,each=nrow(S)))
winner<-rownames(S)[max.col(t(S),ties.method='first')]
eligible<-mx>0 & ties==1L
initial<-data.frame(cell_id=cells$cell_id,stringsAsFactors=FALSE);arms<-list();diagnostics<-list()
for(item in cfg$arms) {
  lambda<-item$lambda;threshold<-lambda*mean(mx)
  calls<-ifelse(eligible & mx>=threshold,winner,'Undecided')
  initial[[item$id]]<-calls
  arms[[item$id]]<-list(arm_id=item$id,library='CM2_glioma_other',cutoff=paste0('cellwise_lambda_',lambda),
    seed_column=item$id,lambda=lambda,threshold=threshold,seed_mechanism='cellwise_positive_normalized_RNA_mean')
  diagnostics[[item$id]]<-data.frame(arm_id=item$id,lambda=lambda,threshold=threshold,
    n_seed=sum(calls!='Undecided'),n_undecided=sum(calls=='Undecided'),n_zero=sum(mx==0),n_tied=sum(ties>1),
    n_classes=length(unique(calls[calls!='Undecided'])))
}
write_atomic(file.path(dest,'initial_calls.csv.gz'),function(p){con<-gzfile(p,'wt');write.csv(initial,con,row.names=FALSE);close(con)})
write_atomic(file.path(dest,'cellwise_scores.rds'),function(p)saveRDS(S,p))
write_atomic(file.path(dest,'cellwise_diagnostics.csv'),function(p)write.csv(data.frame(cell_id=cells$cell_id,max_score=mx,n_ties=ties,positive_unique=eligible,winner=winner),p,row.names=FALSE))
write_atomic(file.path(dest,'seed_arm_counts.csv'),function(p)write.csv(do.call(rbind,diagnostics),p,row.names=FALSE))
write_atomic(file.path(dest,'marker_retention.csv'),function(p)write.csv(do.call(rbind,retention),p,row.names=FALSE))
write_atomic(file.path(dest,'cells.csv'),function(p)write.csv(cells,p,row.names=FALSE))
m<-list(status='score_complete_DL_pending',unit=cfg$sample,sample=cfg$sample,budget=cfg$budget,
  route='cellwise_seed',assay='RNA',scoring_features=nrow(x),n_cells=ncol(x),
  input_signature=cfg$input_signature,protocol_sha256=cfg$protocol_sha256,source_bundle_sha256=cfg$source_bundle_sha256,
  prepare_manifest_sha256=cfg$input$prepare_manifest_sha256,cells_sha256=hash(file.path(dest,'cells.csv')),
  DL_features=pm$features$DL,DL_binary=pm$DL_binary,DL_binary_sha256=pm$DL_binary_sha256,
  DL_input_description='Unchanged native R normalized RNA selected-HVG',
  marker_source_sha256=hash(cfg$marker_file),initial_sha256=hash(file.path(dest,'initial_calls.csv.gz')),
  score_sha256=hash(file.path(dest,'cellwise_scores.rds')),arms=arms,reference_labels_used_for_fit=FALSE,
  no_cluster_or_DEG_used=TRUE,seed_mechanism_replaced=TRUE,mean_cell_max_score=mean(mx),
  elapsed_seconds=as.numeric(difftime(Sys.time(),started,units='secs')),job=Sys.getenv('SLURM_JOB_ID'))
write_atomic(file.path(dest,'score_manifest.json'),function(p)write_json(m,p,auto_unbox=TRUE,pretty=TRUE))
write_atomic(file.path(dest,'SCORE_COMPLETE'),function(p)writeLines(hash(file.path(dest,'score_manifest.json')),p))
cat('A2_SCORE_COMPLETE',cfg$sample,cfg$budget,ncol(x),length(arms),'\n')
