# ScType linear score aggregation, derived from the GPL-3 official function.
# Official source and license remain in vendor_sources; no negative markers
# exist in the shared DG-scRNA panel roster, so gs2=NULL as documented.
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest))
root<-'/fs/scratch/PCON0080/yimin/dgscrna'
base<-file.path(root,'results/hvg_ptc_20260916_v1/paper_claim_validation_20260917')
args<-commandArgs(trailingOnly=TRUE);sample<-args[[1]];mode<-if(length(args)>1L)args[[2]] else 'audit'
prep<-file.path(base,'GBM',sample,'hvg2000')
stopifnot(file.exists(file.path(prep,'PREPARED')))
dest<-file.path(base,'comparators/scType',sample);dir.create(dest,recursive=TRUE,showWarnings=FALSE)
official<-file.path(root,'handoff/paper_claim_validation_20260917/vendor_sources/sctype_score_original.R')
source(official)
libs<-fromJSON(file.path(base,'markers/libraries.json'),simplifyVector=FALSE)
libs<-lapply(libs,function(l)lapply(l,unlist,use.names=FALSE))
obj<-readRDS(file.path(prep,'expression_PCA30.rds'))
genes<-intersect(rownames(obj),unique(unlist(libs,use.names=FALSE)))
obj<-ScaleData(obj,features=genes,do.center=TRUE,do.scale=TRUE,scale.max=10,verbose=FALSE)
Z<-obj[['RNA']]@scale.data
stopifnot(!anyDuplicated(toupper(rownames(Z))))
rownames(Z)<-toupper(rownames(Z))
fast_weights<-function(gs,Z) {
  k<-length(gs)
  stat<-sort(table(unlist(gs)),decreasing=TRUE)
  sensitivity<-setNames(scales::rescale(as.numeric(stat),to=c(0,1),from=c(k,1)),names(stat))
  hit<-lapply(gs,function(g)rownames(Z)[rownames(Z) %in% g])
  keep<-lengths(hit)>0L;hit<-hit[keep]
  i<-rep(seq_along(hit),lengths(hit));j<-match(unlist(hit,use.names=FALSE),rownames(Z))
  w<-unname(sensitivity[rownames(Z)[j]])/rep(sqrt(lengths(hit)),lengths(hit))
  sparseMatrix(i=i,j=j,x=w,dims=c(length(hit),nrow(Z)),dimnames=list(names(hit),rownames(Z)))
}
route_names<-c('PCA30_SNN','PCA30_HDBSCAN_R','UMAP2_SNN','UMAP2_HDBSCAN_R')
is_audit<-startsWith(mode,'audit')
budgets<-if(is_audit)'hvg2000' else c('hvg500','hvg1000','hvg2000','hvg3000','hvg5000','all')
thresholds<-c(default=0.25,permissive=0.125,strict=0.5)
checks<-list();predictions<-data.frame(cell_id=colnames(Z),stringsAsFactors=FALSE);arms<-list()
start<-Sys.time()
for(i in seq_along(libs)) {
  lib<-names(libs)[[i]];W<-fast_weights(libs[[i]],Z)
  if(is_audit) {
    # Smallest prespecified real sample; independently run the unmodified official function.
    original<-sctype_score(scRNAseqData=Z,scaled=TRUE,gs=libs[[i]],gs2=NULL)
    fast<-as.matrix(W %*% Z)
    stopifnot(identical(dimnames(original),dimnames(fast)))
    delta<-max(abs(original-fast))
    stopifnot(is.finite(delta),delta<1e-8)
  }
  for(budget in budgets)for(route in route_names) {
    src<-file.path(base,'GBM',sample,budget,route)
    stopifnot(file.exists(file.path(src,'SCORE_COMPLETE')))
    cl<-read.csv(file.path(src,'clusters.csv'),colClasses='character',stringsAsFactors=FALSE)
    stopifnot(identical(cl$cell_id,colnames(Z)))
    ids<-sort(unique(cl$cluster));idx<-match(cl$cluster,ids)
    A<-sparseMatrix(i=seq_along(idx),j=idx,x=1,dims=c(length(idx),length(ids)))
    sums<-as.matrix(W %*% (Z %*% A));colnames(sums)<-ids
    sizes<-tabulate(idx,nbins=length(ids))
    assign_calls<-function(S,threshold=0.25) {
      winner<-max.col(t(S),ties.method='first')
      values<-S[cbind(winner,seq_len(ncol(S)))]
      result<-rownames(S)[winner];result[values<sizes*threshold]<-'Unknown';result
    }
    if(is_audit) {
      ref<-sapply(seq_along(ids),function(j)rowSums(original[,idx==j,drop=FALSE]))
      rownames(ref)<-rownames(original)
      stopifnot(max(abs(ref-sums))<1e-7)
      for(threshold in thresholds)stopifnot(identical(assign_calls(ref,threshold),assign_calls(sums,threshold)))
      checks[[paste(lib,route,sep='/')]]<-list(cell_score_max_abs_delta=delta,cluster_calls_exact=TRUE)
    }
    for(threshold_name in names(thresholds)) {
    threshold<-thresholds[[threshold_name]];calls<-assign_calls(sums,threshold)
    aid<-paste(budget,route,sprintf('L%02d',i-1L),threshold_name,sep='__')
    predictions[[aid]]<-unname(calls[idx])
    arms[[aid]]<-list(method='scType',library=lib,budget=budget,route=route,threshold=threshold,cutoff=threshold_name,
      n_input_panels=length(libs[[i]]),n_supported_panels=nrow(W),cluster_sha256=digest(file=file.path(src,'clusters.csv'),algo='sha256'))
    }
  }
  cat('SCTYPE_COMPLETE',sample,lib,mode,'\n');flush.console()
}
con<-gzfile(file.path(dest,paste0(mode,'_predictions.csv.gz')),'wt');write.csv(predictions,con,row.names=FALSE);close(con)
m<-list(status='completed',method='scType',sample=sample,mode=mode,arms=arms,checks=checks,n_cells=ncol(Z),
  gene_input='All eligible shared marker genes scaled from full RNA; same per-gene scaling as all-RNA ScaleData',
  no_target_labels_for_fitting=TRUE,negative_markers=NULL,official_source_sha256=digest(file=official,algo='sha256'),
  implementation='Linear sparse aggregation; official per-cell and final cluster call parity audited before cohort use',
  marker_file_sha256=digest(file=file.path(base,'markers/libraries.json'),algo='sha256'),
  elapsed_seconds=as.numeric(difftime(Sys.time(),start,units='secs')),job=Sys.getenv('SLURM_JOB_ID'),
  threshold_grid=thresholds,source_sha256=digest(file=sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]),algo='sha256'),
  predictions_sha256=digest(file=file.path(dest,paste0(mode,'_predictions.csv.gz')),algo='sha256'))
write_json(m,file.path(dest,paste0(mode,'_manifest.json')),pretty=TRUE,auto_unbox=TRUE)
writeLines(digest(file=file.path(dest,paste0(mode,'_manifest.json')),algo='sha256'),file.path(dest,paste0(toupper(mode),'_COMPLETE')))
