#!/usr/bin/env Rscript
# Restore the notebook's two original annotation routes on the archived geometry.
# This is checkpoint-based reconstruction, not a fresh raw-data integration.
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest))
suppressPackageStartupMessages(library(future))
plan(sequential)
root <- '/fs/scratch/PCON0080/yimin/dgscrna'
base <- file.path(root,'results/hvg_ptc_20260916_v1')
arc <- file.path(base,'ptc_recovery/archive/tcr')
task <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
stopifnot(task %in% 0:1)
route <- c('NMT_Thyroid_Seurat_none','TTU_Pubmed_UMAPHDBSCAN_mean')[[task+1L]]
lib <- c('CellMarker_Thyroid','Pubmed_34663816')[[task+1L]]
clname <- c('seurat_clusters','hdbscan.UMAP_clusters')[[task+1L]]
cutoff <- c('none','mean')[[task+1L]]
dest <- file.path(base,'ptc_paper_baseline/replay_selected_routes_marker_union',route)
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
stopifnot(!file.exists(file.path(dest,'SCORE_COMPLETE')))
cat(format(Sys.time()),'Loading historical all-eight-sample checkpoint',route,'\n');flush.console()
checkpoint <- file.path(arc,'rawdata/integrated_data_final_annotation.Rdata')
e <- new.env(); load(checkpoint,envir=e); obj <- e$integrated_data;rm(e);gc()
md <- obj@meta.data
original <- read.delim(file.path(arc,'rawdata/metadata.txt'),check.names=FALSE)
sample_map <- setNames(original$Sample,original$sc_ID)
cells <- paste0(sample_map[md$sample.name],'_',sub('-.*','',rownames(md)))
ref <- read.csv(gzfile(file.path(base,'ptc_paper_baseline/paper_baseline_reference.csv.gz')),
                stringsAsFactors=FALSE,check.names=FALSE)
stopifnot(!anyDuplicated(cells),setequal(cells,ref$cell_id),length(cells)==92404L)
ref <- ref[match(cells,ref$cell_id),,drop=FALSE]
stopifnot(identical(cells,ref$cell_id),all(as.character(md$seurat_clusters)==as.character(ref$seurat_clusters)))
cl <- ref[[if(task==0L)'seurat_clusters' else 'hdbscan_UMAP_clusters']]
ids <- sort(unique(as.integer(cl)))
stopifnot(identical(ids,seq.int(0L,length(ids)-1L)))
obj[[clname]] <- as.character(cl)
Idents(obj) <- factor(as.character(cl),levels=as.character(ids))
DefaultAssay(obj) <- 'integrated'
expr <- obj[['integrated']]@data
stopifnot(nrow(expr)==2000L,ncol(expr)==92404L,identical(colnames(expr),colnames(obj)),
          'CD3D' %in% rownames(expr),all(is.finite(expr)))
writeLines(rownames(expr),file.path(dest,'integration_features.txt'))
write.csv(data.frame(cell_id=cells,original_cell_id=colnames(obj),sample=unname(sample_map[md$sample.name]),
                     group=ref$group,cluster=cl),file.path(dest,'cells.csv'),row.names=FALSE)
# Column-major gene x cell storage is exactly C-order cell x gene for Python.
con <- file(file.path(dest,'DL_archived_CCA2000.float32.bin'),'wb')
for(lo in seq.int(1L,ncol(expr),by=1024L)) {
  hi <- min(ncol(expr),lo+1023L)
  writeBin(as.numeric(expr[,lo:hi,drop=FALSE]),con,size=4,endian='little')
}
close(con)
marker_path <- file.path(arc,'ptc_val/scripts/DGscRNA-Share/data/full_marker_symbol.RDS')
libraries <- readRDS(marker_path)
stopifnot(lib %in% names(libraries))
markers <- libraries[[lib]]
# The recovered later library adds CellMarker_ to the original native prefixes.
# Restore the native vocabulary before fitting; this changes no marker genes/scores.
names(markers) <- sub('^CellMarker_','',names(markers))
write_json(markers,file.path(dest,'marker_panels_original_native_names.json'),pretty=TRUE,auto_unbox=FALSE)
write.csv(data.frame(panel=names(markers),denominator=lengths(markers),
  n_retained=vapply(markers,function(x)length(intersect(x,rownames(expr))),integer(1))),
  file.path(dest,'marker_retention.csv'),row.names=FALSE)
rm(expr);gc()
source_file <- file.path(arc,'ptc_val/scripts/DGscRNA-Share/R/source.R')
# Evaluate the original scoring function without triggering unrelated library loads.
parsed <- parse(source_file)
for(statement in parsed) {
  if(is.call(statement) && identical(statement[[1]],as.name('<-')) &&
     identical(statement[[2]],as.name('density_score'))) eval(statement,envir=.GlobalEnv)
}
stopifnot(exists('density_score'))
legacy_mean <- function(x) log2(Matrix::rowMeans(expm1(x))+1)
deg_file <- file.path(dest,'DEG_marker_union.rds')
if(file.exists(deg_file)) {
  deg <- readRDS(deg_file)
} else {
  cat(format(Sys.time()),'Marker-union FindAllMarkers on unchanged full assay',clname,'clusters',length(ids),'\n');flush.console()
  deg <- FindAllMarkers(obj,assay='integrated',slot='data',test.use='wilcox_limma',
    features=intersect(rownames(obj[['integrated']]),unique(unlist(markers,use.names=FALSE))),
    logfc.threshold=0.25,min.pct=0.1,min.diff.pct=-Inf,only.pos=FALSE,
    max.cells.per.ident=Inf,random.seed=1,min.cells.feature=3,min.cells.group=3,
    pseudocount.use=1,mean.fxn=legacy_mean,fc.name='avg_log2FC',base=2,
    return.thresh=0.01,densify=FALSE,verbose=TRUE)
  stopifnot(nrow(deg)>0L,all(is.finite(deg$avg_log2FC)))
  saveRDS(deg,deg_file)
}
con <- gzfile(file.path(dest,'DEG_marker_union.csv.gz'),'wt');write.csv(deg,con,row.names=FALSE);close(con)
deg_set <- setNames(list(deg),clname)
cat(format(Sys.time()),'Executing archived density_score',lib,cutoff,'\n');flush.console()
obj <- density_score(obj,markers=markers,DEG_markers_set=deg_set,
                     annotation_name=lib,clusterings=clname,cutoffs=cutoff)
col <- paste(lib,clname,cutoff,sep='_')
initial <- as.character(obj@meta.data[[col]])
stopifnot(length(initial)==length(cells),!anyNA(initial))
write.csv(data.frame(cell_id=cells,initial=initial),file.path(dest,'initial_calls.csv'),row.names=FALSE)
write.csv(as.data.frame(table(cluster=cl,initial=initial)),file.path(dest,'cluster_initial_counts.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(dest,'sessionInfo.txt'))
manifest <- list(status='scoring_complete_DL_pending',route=route,library=lib,cutoff=cutoff,cluster_column=clname,
  selected_group=if(task==0L)'NMT' else 'TTU',annotation_fit_scope='all historical 92404 cells; select NMT or TTU after terminal DL',
  geometry='archived all-eight-sample CCA2000 checkpoint; paper S3 partitions; no new integration or cluster fitting',
  sources=list(checkpoint=checkpoint,checkpoint_sha256=digest(file=checkpoint,algo='sha256'),
    source_R_sha256=digest(file=source_file,algo='sha256'),markers_sha256=digest(file=marker_path,algo='sha256'),
    script_sha256=digest(file=file.path(root,'handoff/ptc_paper_baseline_20260916/replay_selected_routes_marker_union.R'),algo='sha256')),
  compatibility='Seurat 5 wilcox_limma and explicit v4 mean/logFC/default thresholds; original density_score executed unmodified',
  optimization='Gene-wise tests restricted to marker union; full assay, normalization, panel denominators and cell contrasts unchanged. Full-gene jobs provide an independent parity check.',scoring_features=2000L,n_cells=length(cells),n_clusters=length(ids),n_initial_undecided=sum(initial=='Undecided'),
  reference_labels_used_for_training=FALSE,job=Sys.getenv('SLURM_JOB_ID'),array_task=task,
  DL_input_sha256=digest(file=file.path(dest,'DL_archived_CCA2000.float32.bin'),algo='sha256'))
write_json(manifest,file.path(dest,'score_manifest.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(digest(file=file.path(dest,'score_manifest.json'),algo='sha256'),file.path(dest,'SCORE_COMPLETE'))
cat(format(Sys.time()),'SCORE_COMPLETE',route,'Undecided',sum(initial=='Undecided'),'\n')
