stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest))
root<-'/fs/scratch/PCON0080/yimin/dgscrna'
base<-file.path(root,'results/hvg_ptc_20260916_v1/r_reference_campaign_20260917')
group<-c('NMT','TTU')[[as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))+1L]]
unit<-paste0('PTC_',group,'_CCA2000');prep<-file.path(base,'PTC_ablation',unit)
original<-file.path(root,'results/hvg_ptc_20260916_v1/ptc_recovery/archive/tcr/ptc_val/scripts/DGscRNA-Share/R/source.R')
for(statement in parse(original))if(is.call(statement)&&identical(statement[[1]],as.name('<-'))&&identical(statement[[2]],as.name('density_score')))eval(statement,envir=.GlobalEnv)
libs<-fromJSON(file.path(base,'markers/PTC_original17.json'),simplifyVector=FALSE)
libs<-lapply(libs,function(l)lapply(l,unlist,use.names=FALSE))
checks<-list();diagnostics<-list()
for(route in c('PCA30_SNN','PCA30_HDBSCAN_R','UMAP2_SNN','UMAP2_HDBSCAN_R')) {
  source<-file.path(prep,route)
  cl<-read.csv(file.path(source,'clusters.csv'),stringsAsFactors=FALSE,colClasses='character')
  ids<-sort(unique(cl$cluster));index<-setNames(seq_along(ids)-1L,ids)
  deg<-readRDS(file.path(source,'DEG.rds'));deg$cluster<-unname(index[as.character(deg$cluster)])
  initial<-read.csv(gzfile(file.path(source,'initial_calls.csv.gz')),stringsAsFactors=FALSE,check.names=FALSE,colClasses='character')
  stopifnot(identical(initial$cell_id,cl$cell_id))
  dummy<-matrix(1,nrow=2,ncol=length(ids),dimnames=list(c('AUDIT_A','AUDIT_B'),paste0('cluster_',ids)))
  obj<-CreateSeuratObject(dummy,min.cells=0,min.features=0)
  obj$reference_clusters<-as.character(seq_along(ids)-1L)
  for(i in seq_along(libs)) {
    lib<-names(libs)[[i]]
    obj<-density_score(obj,markers=libs[[i]],DEG_markers_set=list(reference_clusters=deg),annotation_name=lib,
                        clusterings='reference_clusters',cutoffs=c('none','mean','0.5'))
    for(cut in c('none','mean','0.5')) {
      aid<-sprintf('L%02d_%s',i-1L,if(cut=='0.5')'p050' else cut)
      metadata_name<-gsub('[ /]','.',paste(lib,'reference_clusters',cut,sep='_'))
      stopifnot(metadata_name %in% colnames(obj@meta.data))
      calls<-obj@meta.data[[metadata_name]]
      got<-as.character(calls[unname(index[cl$cluster])+1L])
      if(!identical(got,as.character(initial[[aid]]))) {
        cat('MISMATCH',unit,route,lib,cut,'length',length(got),'expected',length(initial[[aid]]),'\n')
        print(head(data.frame(replayed=got,saved=initial[[aid]])[got!=initial[[aid]],],10))
        stop('PTC original-density parity failed')
      }
      checks[[paste(route,aid,sep='/')]]<-list(exact=TRUE,n_cells=nrow(cl),library=lib,cutoff=cut)
    }
  }
  lib<-if(group=='NMT')'CellMarker_Thyroid' else 'Pubmed_34663816'
  tnames<-if(group=='NMT')'cancer+Thyroid+Thyroid+T cell' else c('NCOMMREFF+T Cells','NCOMMREFF+Treg Cells')
  S<-readRDS(file.path(source,'density_scores.rds'))[[lib]]
  stopifnot(all(tnames %in% rownames(S)))
  maximum<-apply(S,2,max);threshold<-mean(maximum)
  wins<-vapply(colnames(S),function(k) {
    winners<-rownames(S)[S[,k]==maximum[[k]]]
    length(winners)==1L && winners %in% tnames
  },logical(1))
  diagnostics[[route]]<-list(library=lib,T_panels=tnames,mean_cutoff=threshold,
    T_winning_clusters_before_cutoff=sum(wins),T_winning_clusters_after_mean=sum(wins & maximum>=threshold),
    T_winning_cells_before_cutoff=sum(cl$cluster %in% names(wins)[wins]),
    T_winning_cells_after_mean=sum(cl$cluster %in% names(wins)[wins & maximum>=threshold]),
    T_panel_maximum_scores=apply(S[tnames,,drop=FALSE],1,max),
    cluster_T_panel_scores=S[tnames,,drop=FALSE])
  cat('PARITY',unit,route,length(libs)*3L,'\n');flush.console()
}
dest<-file.path(base,'verification/ptc_density');dir.create(dest,recursive=TRUE,showWarnings=FALSE)
write_json(list(status='passed',unit=unit,job=Sys.getenv('SLURM_JOB_ID'),checks=checks,diagnostics=diagnostics,
   original_function_source_sha256=digest(file=original,algo='sha256'),
   script_sha256=digest(file=file.path(root,'handoff/r_reference_campaign_20260917/audit_ptc_density.R'),algo='sha256'),
   meaning='All 204 initial annotation arms exactly replay the original density function; diagnostic scores distinguish missing markers from cutoff/competition.'),
   file.path(dest,paste0(unit,'.json')),auto_unbox=TRUE,pretty=TRUE)
cat(toJSON(diagnostics,auto_unbox=TRUE,pretty=TRUE),'\n')
