#!/usr/bin/env Rscript
# Reconstruct the historical density rule with explicit Seurat-v4 statistics.
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest))
suppressPackageStartupMessages(library(future))
plan(sequential);options(Seurat.object.assay.version='v3')
args <- commandArgs(trailingOnly=TRUE)
mode <- args[[1]]
task_id <- if(length(args)>1) as.integer(args[[2]]) else as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
root <- '/fs/scratch/PCON0080/yimin/dgscrna'
base <- file.path(root,'results/hvg_ptc_20260916_v1/ptc_experiments')
marker_path <- file.path(root,'results/hvg_ptc_20260916_v1/ptc_recovery/archive/tcr/ptc_val/scripts/DGscRNA-Share/data/full_marker_symbol.RDS')
libraries <- readRDS(marker_path)
stopifnot(length(libraries)==17L)
legacy_mean <- function(x) log2(Matrix::rowMeans(expm1(x))+1)
source_path <- Sys.getenv('PTC_SCORE_SOURCE')
stopifnot(nzchar(source_path))
source_sha <- digest(file=source_path,algo='sha256')

write_gz <- function(x,path,row.names=FALSE) {
  con <- gzfile(paste0(path,'.part'),'wt');write.csv(x,con,row.names=row.names);close(con)
  stopifnot(file.rename(paste0(path,'.part'),path))
}

score_one <- function(obj,pdir,assay,selected_libraries,cutoffs,context,parity=FALSE) {
  stopifnot(file.exists(file.path(pdir,'CLUSTER_COMPLETE')))
  dest <- file.path(pdir,paste0('score_',assay));dir.create(dest,showWarnings=FALSE)
  if(file.exists(file.path(dest,'SCORE_COMPLETE'))) {
    previous <- fromJSON(file.path(dest,'score_manifest.json'),simplifyVector=FALSE)
    stopifnot(identical(previous$source_sha256,source_sha))
    return(invisible(NULL))
  }
  started <- Sys.time()
  cl <- read.csv(file.path(pdir,'clusters.csv'),stringsAsFactors=FALSE)
  stopifnot(!anyDuplicated(cl$cell_id),setequal(cl$cell_id,colnames(obj)))
  cl <- cl[match(colnames(obj),cl$cell_id),,drop=FALSE]
  stopifnot(identical(cl$cell_id,colnames(obj)))
  ids <- sort(unique(as.character(cl$cluster)))
  Idents(obj) <- factor(as.character(cl$cluster),levels=ids)
  DefaultAssay(obj) <- assay
  # Each Wilcoxon test and its fold change is gene-wise. Only genes that can enter
  # these panels need testing; normalization and the full assay universe are unchanged.
  marker_union <- unique(unlist(libraries[selected_libraries],use.names=FALSE))
  features <- intersect(rownames(obj[[assay]]),marker_union)
  stopifnot(length(features)>0)
  warnings_seen <- character();skipped <- list()
  deg_path <- file.path(dest,'DEG_marker_union.rds')
  if(file.exists(deg_path)) {
    cached <- readRDS(deg_path)
    stopifnot(cached$cluster_sha256==digest(file=file.path(pdir,'clusters.csv'),algo='sha256'),
              cached$source_sha256==source_sha)
    deg <- cached$deg;skipped <- cached$skipped;warnings_seen <- cached$warnings
  } else {
    degs <- list()
    if(length(ids)>1L) for(id in ids) {
      n_in <- sum(Idents(obj)==id);n_out <- ncol(obj)-n_in
      if(min(n_in,n_out)<3L) {
        skipped[[id]] <- paste0('structural: fewer than 3 cells in a contrast (',n_in,',',n_out,')')
        next
      }
      cat(format(Sys.time(),tz='UTC'),'DEG',pdir,assay,id,n_in,'\n');flush.console()
      d <- withCallingHandlers(FindMarkers(obj,assay=assay,ident.1=id,ident.2=NULL,
        features=features,slot='data',test.use='wilcox_limma',logfc.threshold=0.25,
        min.pct=0.1,min.diff.pct=-Inf,only.pos=FALSE,max.cells.per.ident=Inf,random.seed=1,
        min.cells.feature=3,min.cells.group=3,pseudocount.use=1,mean.fxn=legacy_mean,
        fc.name='avg_log2FC',base=2,densify=FALSE,verbose=FALSE),warning=function(w) {
          warnings_seen <<- c(warnings_seen,paste(id,conditionMessage(w)))
      })
      if(nrow(d)>0L) {
        stopifnot(all(c('p_val','avg_log2FC') %in% colnames(d)),
                  all(is.finite(d$avg_log2FC)),all(is.finite(d$p_val)))
        d$gene <- rownames(d);d$cluster <- id
        d <- d[d$p_val<0.01,,drop=FALSE]
        if(nrow(d)>0L) degs[[id]] <- d
      }
    }
    deg <- if(length(degs)) do.call(rbind,degs) else data.frame(p_val=numeric(),
            avg_log2FC=numeric(),pct.1=numeric(),pct.2=numeric(),p_val_adj=numeric(),
            gene=character(),cluster=character())
    saveRDS(list(deg=deg,skipped=skipped,warnings=warnings_seen,
      cluster_sha256=digest(file=file.path(pdir,'clusters.csv'),algo='sha256'),
      source_sha256=source_sha),deg_path)
  }
  write_gz(deg,file.path(dest,'DEG_marker_union.csv.gz'))
  if(parity) {
    cat('Independent full-gene FindAllMarkers parity check\n');flush.console()
    full <- FindAllMarkers(obj,assay=assay,slot='data',test.use='wilcox_limma',
      logfc.threshold=0.25,min.pct=0.1,min.diff.pct=-Inf,only.pos=FALSE,
      max.cells.per.ident=Inf,random.seed=1,min.cells.feature=3,min.cells.group=3,
      pseudocount.use=1,mean.fxn=legacy_mean,fc.name='avg_log2FC',base=2,
      return.thresh=0.01,densify=FALSE,verbose=FALSE)
    write_gz(full,file.path(dest,'parity_full_gene_FindAllMarkers.csv.gz'))
    projected <- full[full$gene %in% features,,drop=FALSE]
    key <- function(d) paste(d$cluster,d$gene,sep='||')
    stopifnot(setequal(key(projected),key(deg)))
    projected <- projected[match(key(deg),key(projected)),,drop=FALSE]
    for(n in c('p_val','avg_log2FC','pct.1','pct.2','p_val_adj'))
      stopifnot(isTRUE(all.equal(projected[[n]],deg[[n]],tolerance=1e-12,check.attributes=FALSE)))
    write_json(list(status='passed',n_significant_panel_genes=nrow(deg),
      full_gene_count=nrow(obj[[assay]]),tested_union_count=length(features),
      criterion='Full-gene FindAllMarkers projected onto panel union equals per-cluster union tests, all five statistics to 1e-12'),
      file.path(dest,'full_gene_parity.json'),pretty=TRUE,auto_unbox=TRUE)
  }
  calls <- list();scores <- list();arms <- list();retention <- list()
  cell_seed <- data.frame(cell_id=colnames(obj),stringsAsFactors=FALSE)
  for(lib in selected_libraries) {
    panels <- libraries[[lib]]
    stopifnot(all(lengths(panels)>0L),!anyDuplicated(names(panels)))
    S <- matrix(0,nrow=length(panels),ncol=length(ids),dimnames=list(names(panels),ids))
    for(id in ids) {
      ds <- deg[as.character(deg$cluster)==id & deg$avg_log2FC>1,,drop=FALSE]
      for(panel in names(panels)) {
        hits <- intersect(panels[[panel]],ds$gene)
        S[panel,id] <- sum(ds$avg_log2FC[ds$gene %in% hits])/length(panels[[panel]])
        if(length(panels[[panel]])<=1L) S[panel,id] <- S[panel,id]*0.8
      }
    }
    scores[[lib]] <- S
    retention[[lib]] <- data.frame(library=lib,panel=names(panels),
      full_denominator=lengths(panels),unique_symbols=sapply(panels,function(x)length(unique(x))),
      retained_unique=sapply(panels,function(x)length(intersect(x,rownames(obj[[assay]])))))
    mx <- apply(S,2,max)
    for(cutoff in cutoffs) {
      cluster_calls <- vapply(ids,function(id) {
        winners <- rownames(S)[S[,id]==mx[[id]]]
        answer <- if(length(winners)==1L) winners else 'Undecided'
        threshold <- if(cutoff=='mean') mean(mx) else if(cutoff=='0.5') 0.5 else -Inf
        if(mx[[id]]<threshold) answer <- 'Undecided'
        answer
      },character(1))
      aid <- sprintf('L%02d_%s',match(lib,names(libraries))-1L,if(cutoff=='0.5')'p050' else cutoff)
      cell_seed[[aid]] <- unname(cluster_calls[as.character(cl$cluster)])
      stopifnot(!anyNA(cell_seed[[aid]]))
      calls[[aid]] <- data.frame(arm_id=aid,library=lib,cutoff=cutoff,cluster=ids,
                  max_score=as.numeric(mx),initial=unname(cluster_calls))
      arms[[aid]] <- list(arm_id=aid,library=lib,cutoff=cutoff,seed_column=aid,
                        scoring_assay=assay,DL_input=if(assay=='integrated')'CCA2000' else 'RNA2000')
    }
  }
  saveRDS(scores,file.path(dest,'panel_scores.rds'))
  write_gz(cell_seed,file.path(dest,'initial_calls.csv.gz'))
  write_gz(do.call(rbind,calls),file.path(dest,'cluster_calls.csv.gz'))
  write_gz(do.call(rbind,retention),file.path(dest,'marker_retention.csv.gz'))
  m <- list(status='completed',context=context,assay=assay,n_cells=ncol(obj),n_clusters=length(ids),
      full_scoring_gene_universe=nrow(obj[[assay]]),tested_marker_union=length(features),
      DEG_optimization='Only marker-union genes need gene-wise tests; unchanged full RNA normalization and panel denominators. Independently verified by full-gene pilot.',
      DE_rule=list(test='wilcox_limma',logfc_threshold=0.25,min_pct=0.1,
        raw_p_return_threshold=0.01,scoring_log2FC_threshold=1,
        mean_fxn='log2(rowMeans(expm1(x))+1)',singleton_penalty=0.8),
      arms=arms,source_sha256=source_sha,marker_sha256=digest(file=marker_path,algo='sha256'),
      cluster_sha256=digest(file=file.path(pdir,'clusters.csv'),algo='sha256'),
      initial_sha256=digest(file=file.path(dest,'initial_calls.csv.gz'),algo='sha256'),
      observed_cluster_ID_mapping=as.list(setNames(ids,seq_along(ids)-1L)),
      cluster_index_adapter='Iterate actual observed IDs, equivalent to explicit contiguous 0-based relabeling; does not drop highest label when historical loop assumptions fail',
      structural_contrasts=skipped,warnings=unique(warnings_seen),
      R=R.version.string,Seurat=as.character(packageVersion('Seurat')),
      limma=if(requireNamespace('limma',quietly=TRUE))as.character(packageVersion('limma')) else NA_character_,
      job=Sys.getenv('SLURM_JOB_ID'),elapsed_seconds=as.numeric(difftime(Sys.time(),started,units='secs')))
  write_json(m,file.path(dest,'score_manifest.json'),pretty=TRUE,auto_unbox=TRUE,na='null')
  writeLines(digest(file=file.path(dest,'score_manifest.json'),algo='sha256'),file.path(dest,'SCORE_COMPLETE'))
  invisible(NULL)
}

if(mode %in% c('single','pilot')) {
  tasks <- fromJSON(file.path(base,'protocol/single_sample_geometries.json'),simplifyVector=FALSE)
  t <- tasks[[task_id+1L]]
  obj <- readRDS(file.path(base,'samples',t$sample,'RNA_normalized_VST.rds'))
  gd <- file.path(base,'single_sample',t$sample,t$geometry_id)
  stopifnot(file.exists(file.path(gd,'GEOMETRY_COMPLETE')))
  libs <- if(isTRUE(t$full_roster)) names(libraries) else 'CellMarker_AllTissues'
  cutoffs <- if(isTRUE(t$full_roster)) c('none','mean','0.5') else 'mean'
  for(method in t$clusterers) {
    if(mode=='pilot' && method!='HDBSCAN') next
    score_one(obj,file.path(gd,method),'RNA',libs,cutoffs,t,parity=mode=='pilot')
  }
} else if(mode=='pooled') {
  group <- c('MTN','TUT')[[task_id %/% 3L+1L]]
  correction <- c('NONE','CCA','HARMONY')[[task_id %% 3L+1L]]
  gd <- file.path(base,'pooled',group,correction)
  stopifnot(file.exists(file.path(gd,'GEOMETRY_COMPLETE')))
  obj <- readRDS(file.path(base,'prepared',group,'RNA_normalized.rds'))
  for(space in c('PCA30','UMAP2')) for(method in c('SNN','HDBSCAN_R')) {
    context <- list(group=group,correction=correction,space=space,clusterer=method,seed=42)
    score_one(obj,file.path(gd,paste0(space,'_',method)),'RNA',names(libraries),c('none','mean','0.5'),context)
  }
  rm(obj);gc()
  if(correction=='CCA') {
    obj <- readRDS(file.path(base,'prepared',group,'CCA_PCA30.rds'))
    for(space in c('PCA30','UMAP2')) for(method in c('SNN','HDBSCAN_R')) {
      context <- list(group=group,correction=correction,space=space,clusterer=method,seed=42,
                     family='legacy_integrated_scoring_and_DL')
      score_one(obj,file.path(gd,paste0(space,'_',method)),'integrated',names(libraries),c('none','mean','0.5'),context)
    }
  }
} else stop('Unknown execution mode')
cat('All requested R scoring partitions completed\n')
