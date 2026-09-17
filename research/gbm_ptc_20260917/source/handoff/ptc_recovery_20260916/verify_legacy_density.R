#!/usr/bin/env Rscript
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest))
root <- '/fs/scratch/PCON0080/yimin/dgscrna'
base <- file.path(root,'results/hvg_ptc_20260916_v1/ptc_experiments')
archive <- file.path(root,'results/hvg_ptc_20260916_v1/ptc_recovery/archive/tcr/ptc_val/scripts/DGscRNA-Share')
source_path <- file.path(archive,'R/source.R')
exprs <- parse(source_path)
target <- which(vapply(exprs,function(e) is.call(e) && identical(e[[1]],as.name('<-')) &&
                         identical(e[[2]],as.name('density_score')),logical(1)))
stopifnot(length(target)==1L)
eval(exprs[[target]])
tasks <- fromJSON(file.path(base,'protocol/single_sample_geometries.json'),simplifyVector=FALSE)
t <- tasks[[1]]
pdir <- file.path(base,'single_sample',t$sample,t$geometry_id,'HDBSCAN')
scored <- file.path(pdir,'score_RNA')
stopifnot(file.exists(file.path(scored,'full_gene_parity.json')))
obj <- readRDS(file.path(base,'samples',t$sample,'RNA_normalized_VST.rds'))
cl <- read.csv(file.path(pdir,'clusters.csv'),stringsAsFactors=FALSE)
cl <- cl[match(colnames(obj),cl$cell_id),,drop=FALSE]
ids <- sort(unique(as.character(cl$cluster)))
map <- setNames(seq_along(ids)-1L,ids)
obj$louvain <- factor(unname(map[as.character(cl$cluster)]))
deg <- readRDS(file.path(scored,'DEG_marker_union.rds'))$deg
deg$cluster <- unname(map[as.character(deg$cluster)])
libraries <- readRDS(file.path(archive,'data/full_marker_symbol.RDS'))
manifest <- fromJSON(file.path(scored,'score_manifest.json'),simplifyVector=FALSE)
initial <- read.csv(gzfile(file.path(scored,'initial_calls.csv.gz')),check.names=FALSE,
                    stringsAsFactors=FALSE,na.strings=NULL)
initial <- initial[match(colnames(obj),initial$cell_id),,drop=FALSE]
for(arm in manifest$arms) {
  result <- density_score(obj,markers=libraries[[arm$library]],DEG_markers_set=list(louvain=deg),
                          annotation_name='Legacy',clusterings='louvain',cutoffs=arm$cutoff)
  column <- paste0('Legacy_louvain_',arm$cutoff)
  stopifnot(identical(as.character(result[[column]][,1]),initial[[arm$arm_id]]))
}
dest <- file.path(base,'verification/legacy_density');dir.create(dest,recursive=TRUE,showWarnings=FALSE)
write_json(list(status='passed',original_source_sha256=digest(file=source_path,algo='sha256'),
   scorer_sha256=manifest$source_sha256,n_cells=ncol(obj),n_arms=length(manifest$arms),
   check='Original archived density_score function, evaluated unchanged, agrees per cell after explicit contiguous cluster-index adapter'),
   file.path(dest,'verification.json'),pretty=TRUE,auto_unbox=TRUE)
writeLines(digest(file=file.path(dest,'verification.json'),algo='sha256'),file.path(dest,'COMPLETE'))
cat('Archived R density function agrees exactly\n')
