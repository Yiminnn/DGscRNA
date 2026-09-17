#!/usr/bin/env Rscript
# Read the original saved R workspace in an isolated environment and export metadata.
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
stopifnot(requireNamespace('SeuratObject',quietly=TRUE),requireNamespace('jsonlite',quietly=TRUE))
root <- '/fs/scratch/PCON0080/yimin/dgscrna/results/hvg_ptc_20260916_v1/ptc_recovery'
stopifnot(file.exists(file.path(root,'ARCHIVE_STAGED')))
source_path <- file.path(root,'archive/tcr/rawdata/integrated_data_final_annotation.Rdata')
dest <- file.path(root,'inventory_workspace')
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
historical <- new.env(parent=emptyenv())
cat(format(Sys.time(),tz='UTC'),'Loading saved R workspace\n');flush.console()
loaded <- load(source_path,envir=historical)
writeLines(loaded,file.path(dest,'object_names.txt'))
writeLines(capture.output(sessionInfo()),file.path(dest,'sessionInfo.txt'))
inventory <- list()
for (i in seq_along(loaded)) {
  name <- loaded[[i]]
  obj <- get(name,envir=historical,inherits=FALSE)
  key <- sprintf('object_%02d',i)
  info <- list(name=name,key=key,class=class(obj),dim=dim(obj),
               size_bytes=as.numeric(object.size(obj)))
  if (inherits(obj,'Seurat')) {
    meta <- methods::slot(obj,'meta.data')
    write.csv(meta,file.path(dest,paste0(key,'.meta.csv')))
    info$active_assay <- methods::slot(obj,'active.assay')
    info$assays <- lapply(methods::slot(obj,'assays'),function(a) list(class=class(a),dim=dim(a),slots=methods::slotNames(a)))
    info$reductions <- lapply(methods::slot(obj,'reductions'),function(x) dim(methods::slot(x,'cell.embeddings')))
    info$seurat_version <- as.character(methods::slot(obj,'version'))
    commands <- methods::slot(obj,'commands')
    info$commands <- lapply(commands,function(x) list(name=methods::slot(x,'name'),
      assay=methods::slot(x,'assay.used'),call=methods::slot(x,'call.string'),
      time=as.character(methods::slot(x,'time.stamp'))))
    writeLines(capture.output(str(commands,max.level=4,list.len=60,give.attr=FALSE)),
      file.path(dest,paste0(key,'.commands.txt')))
    for (assay_name in names(methods::slot(obj,'assays'))) {
      assay <- methods::slot(obj,'assays')[[assay_name]]
      if ('var.features' %in% methods::slotNames(assay))
        writeLines(as.character(methods::slot(assay,'var.features')),
          file.path(dest,paste0(key,'.',assay_name,'.variable_features.txt')))
    }
    info$meta_columns <- colnames(meta)
    info$meta_low_cardinality <- lapply(meta[vapply(meta,function(x) length(unique(x))<=200,logical(1))],
                                        function(x) as.list(table(x,useNA='ifany')))
    write.csv(data.frame(gene=rownames(obj)),file.path(dest,paste0(key,'.genes.csv')),row.names=FALSE)
    writeLines(capture.output(str(obj,max.level=2,list.len=20,give.attr=FALSE)),file.path(dest,paste0(key,'.structure.txt')))
  } else if (is.data.frame(obj)) {
    info$columns <- colnames(obj)
    write.csv(obj,file.path(dest,paste0(key,'.table.csv')))
  } else if (is.function(obj)) {
    # Saved workspaces can retain the exact historical helper bodies.
    writeLines(deparse(obj,width.cutoff=120),file.path(dest,paste0(key,'.function.R')))
    info$formals <- names(formals(obj))
  } else if (is.list(obj)) {
    info$length <- length(obj)
    info$names <- names(obj)
    info$children <- lapply(obj,function(x) list(class=class(x),dim=dim(x),length=length(x)))
    writeLines(capture.output(str(obj,max.level=2,list.len=30,give.attr=FALSE)),file.path(dest,paste0(key,'.structure.txt')))
  } else {
    writeLines(capture.output(str(obj,max.level=2,list.len=30,give.attr=FALSE)),file.path(dest,paste0(key,'.structure.txt')))
  }
  inventory[[name]] <- info
  jsonlite::write_json(list(source=source_path,job=Sys.getenv('SLURM_JOB_ID'),objects=inventory),
                      file.path(dest,'objects.json'),pretty=TRUE,auto_unbox=TRUE,na='null')
  cat('Inventoried',name,paste(class(obj),collapse='/'),info$size_bytes,'bytes\n');flush.console()
  rm(obj);gc()
}
writeLines(format(Sys.time(),tz='UTC'),file.path(dest,'COMPLETE'))
