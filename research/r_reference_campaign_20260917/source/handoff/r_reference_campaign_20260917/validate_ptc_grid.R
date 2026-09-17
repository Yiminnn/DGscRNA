stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest))
root<-'/fs/scratch/PCON0080/yimin/dgscrna'
base<-file.path(root,'results/hvg_ptc_20260916_v1')
out<-file.path(base,'r_reference_campaign_20260917')
libs<-readRDS(file.path(base,'ptc_recovery/archive/tcr/ptc_val/scripts/DGscRNA-Share/data/full_marker_symbol.RDS'))
for(lib in names(libs)) names(libs[[lib]])<-sub('^CellMarker_','',names(libs[[lib]]))
write_json(libs,file.path(out,'markers/PTC_original17.json'),auto_unbox=FALSE,pretty=TRUE)
checks<-list()
for(spec in list(c('seurat_clusters','CellMarker_Thyroid','none','NMT_Thyroid_Seurat_none'),
                c('hdbscan.UMAP_clusters','Pubmed_34663816','mean','TTU_Pubmed_UMAPHDBSCAN_mean'))) {
  dest<-file.path(out,'PTC_archived_CCA2000',spec[[1]])
  if(!file.exists(file.path(dest,'SCORE_COMPLETE')))next
  m<-fromJSON(file.path(dest,'score_manifest.json'),simplifyVector=FALSE)
  wanted<-Filter(function(a)a$library==spec[[2]] && a$cutoff==spec[[3]],m$arms)
  stopifnot(length(wanted)==1L)
  calls<-read.csv(gzfile(file.path(dest,'initial_calls.csv.gz')),check.names=FALSE,stringsAsFactors=FALSE)
  original<-read.csv(file.path(base,'ptc_paper_baseline/replay_selected_routes_full_parallel',spec[[4]],'initial_calls.csv'),stringsAsFactors=FALSE)
  stopifnot(identical(calls$cell_id,original$cell_id),identical(calls[[wanted[[1]]$seed_column]],original$initial))
  checks[[spec[[4]]]]<-list(initial_calls_exact=TRUE,n_cells=nrow(calls),arm=wanted[[1]]$arm_id,
    score_manifest_sha256=digest(file=file.path(dest,'score_manifest.json'),algo='sha256'))
}
write_json(list(checks=checks,job=Sys.getenv('SLURM_JOB_ID'),libraries=length(libs)),file.path(out,'ptc_selected_initial_parity.json'),auto_unbox=TRUE,pretty=TRUE)
cat('PTC marker export and available selected-route parity passed\n')
