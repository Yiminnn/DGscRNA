#!/usr/bin/env Rscript
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest))
suppressPackageStartupMessages(library(future))
plan(multicore,workers=as.integer(Sys.getenv('SLURM_CPUS_PER_TASK','4')))
options(future.globals.maxSize=16*1024^3)
root <- '/fs/scratch/PCON0080/yimin/dgscrna'
base <- file.path(root,'results/hvg_ptc_20260916_v1')
arc <- file.path(base,'ptc_recovery/archive/tcr')
out <- file.path(base,'r_reference_campaign_20260917/PTC_archived_CCA2000')
task <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
cluster_columns <- c('seurat_clusters','seurat.UMAP_clusters','hdbscan_clusters','hdbscan.UMAP_clusters')
clname <- cluster_columns[[task+1L]]
dest <- file.path(out,clname);dir.create(dest,recursive=TRUE,showWarnings=FALSE)
source_sha <- digest(file=file.path(root,'handoff/r_reference_campaign_20260917/ptc_grid.R'),algo='sha256')
if(file.exists(file.path(dest,'SCORE_COMPLETE'))) {
  old <- fromJSON(file.path(dest,'score_manifest.json'))
  stopifnot(old$source_sha256==source_sha)
  quit(status=0)
}
start <- Sys.time()
e <- new.env();load(file.path(arc,'rawdata/integrated_data_final_annotation.Rdata'),envir=e)
obj <- e$integrated_data;rm(e);gc()
previous <- file.path(base,'ptc_paper_baseline/replay_selected_routes_full_parallel/NMT_Thyroid_Seurat_none')
cells <- read.csv(file.path(previous,'cells.csv'),stringsAsFactors=FALSE,check.names=FALSE)
stopifnot(identical(cells$original_cell_id,colnames(obj)))
original <- read.csv(file.path(arc,'rawdata/data_with_validation.csv'),stringsAsFactors=FALSE,check.names=FALSE)
sample_metadata <- read.delim(file.path(arc,'rawdata/metadata.txt'),check.names=FALSE)
sample_map <- setNames(sample_metadata$Sample,sample_metadata$sc_ID)
original_ids <- paste0(sample_map[as.character(original$sample.name)],'_',sub('-.*','',original[[1]]))
stopifnot(!anyDuplicated(original_ids),setequal(original_ids,cells$cell_id),clname %in% colnames(original))
original <- original[match(cells$cell_id,original_ids),,drop=FALSE]
cl <- as.character(original[[clname]])
obj[[clname]] <- cl
stopifnot(identical(sort(unique(as.integer(cl))),seq.int(0L,length(unique(cl))-1L)))
stopifnot(!anyNA(cl))
cells$cluster <- cl
write.csv(cells,file.path(dest,'cells.csv'),row.names=FALSE)
Idents(obj) <- factor(cl,levels=as.character(sort(unique(as.integer(cl)))))
DefaultAssay(obj) <- 'integrated'
libraries <- readRDS(file.path(arc,'ptc_val/scripts/DGscRNA-Share/data/full_marker_symbol.RDS'))
stopifnot(length(libraries)==17L)
for(lib in names(libraries)) names(libraries[[lib]]) <- sub('^CellMarker_','',names(libraries[[lib]]))
parsed <- parse(file.path(arc,'ptc_val/scripts/DGscRNA-Share/R/source.R'))
for(statement in parsed) if(is.call(statement) && identical(statement[[1]],as.name('<-')) &&
    identical(statement[[2]],as.name('density_score'))) eval(statement,envir=.GlobalEnv)
legacy_mean <- function(x) log2(Matrix::rowMeans(expm1(x))+1)
cached <- switch(clname,
  seurat_clusters=file.path(previous,'DEG_full_integrated2000.rds'),
  hdbscan.UMAP_clusters=file.path(base,'ptc_paper_baseline/replay_selected_routes_full_parallel/TTU_Pubmed_UMAPHDBSCAN_mean/DEG_full_integrated2000.rds'),
  file.path(dest,'DEG_full_integrated2000.rds'))
if(file.exists(cached)) {
  deg <- readRDS(cached)
} else {
  deg <- FindAllMarkers(obj,assay='integrated',slot='data',test.use='wilcox_limma',
    logfc.threshold=0.25,min.pct=0.1,min.diff.pct=-Inf,only.pos=FALSE,
    max.cells.per.ident=Inf,random.seed=1,min.cells.feature=3,min.cells.group=3,
    pseudocount.use=1,mean.fxn=legacy_mean,fc.name='avg_log2FC',base=2,
    return.thresh=0.01,densify=FALSE,verbose=TRUE)
  stopifnot(nrow(deg)>0L,all(is.finite(deg$avg_log2FC)))
  saveRDS(deg,cached)
}
degs <- setNames(list(deg),clname)
initial <- data.frame(cell_id=cells$cell_id,stringsAsFactors=FALSE)
arms <- list();retention <- list()
for(i in seq_along(libraries)) {
  lib <- names(libraries)[[i]];panels <- libraries[[i]]
  obj <- density_score(obj,markers=panels,DEG_markers_set=degs,annotation_name=lib,
                       clusterings=clname,cutoffs=c('none','mean','0.5'))
  retention[[lib]] <- data.frame(library=lib,panel=names(panels),denominator=lengths(panels),
    retained=vapply(panels,function(g)length(intersect(g,rownames(obj[['integrated']]))),integer(1)))
  for(cut in c('none','mean','0.5')) {
    aid <- sprintf('L%02d_%s',i-1L,if(cut=='0.5')'p050' else cut)
    metadata_column <- gsub('[ /]','.',paste(lib,clname,cut,sep='_'))
    stopifnot(metadata_column %in% colnames(obj@meta.data))
    initial[[aid]] <- as.character(obj@meta.data[[metadata_column]])
    stopifnot(!anyNA(initial[[aid]]))
    arms[[aid]] <- list(arm_id=aid,library=lib,cutoff=cut,seed_column=aid)
  }
}
con <- gzfile(file.path(dest,'initial_calls.csv.gz'),'wt');write.csv(initial,con,row.names=FALSE);close(con)
write.csv(do.call(rbind,retention),file.path(dest,'marker_retention.csv'),row.names=FALSE)
binary <- file.path(previous,'DL_archived_CCA2000.float32.bin')
m <- list(status='score_complete_DL_pending',dataset='PTC',fit_scope='all original 92404 cells',
  report_groups=c('NMT','TTU'),geometry='archived all-eight-sample CCA/PCA30/UMAP2',
  cluster_column=clname,source_sha256=source_sha,n_cells=ncol(obj),scoring_features=2000L,
  DL_features=2000L,DL_binary=binary,DL_binary_sha256=digest(file=binary,algo='sha256'),
  initial_sha256=digest(file=file.path(dest,'initial_calls.csv.gz'),algo='sha256'),
  cell_order_sha256=digest(file=file.path(dest,'cells.csv'),algo='sha256'),
  DEG_file=cached,DEG_sha256=digest(file=cached,algo='sha256'),
  original_density_function=TRUE,noise_zero_scored_as_cluster=TRUE,arms=arms,
  annotation_labels_used_for_fit=FALSE,job=Sys.getenv('SLURM_JOB_ID'),
  elapsed_seconds=as.numeric(difftime(Sys.time(),start,units='secs')))
write_json(m,file.path(dest,'score_manifest.json'),auto_unbox=TRUE,pretty=TRUE)
writeLines(capture.output(sessionInfo()),file.path(dest,'sessionInfo.txt'))
writeLines(digest(file=file.path(dest,'score_manifest.json'),algo='sha256'),file.path(dest,'SCORE_COMPLETE'))
cat('COMPLETE',clname,length(arms),'annotation arms\n')
