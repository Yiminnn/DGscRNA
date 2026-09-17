#!/usr/bin/env Rscript
# Inventory historical marker libraries and R objects. No preprocessing or fitting.
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
root <- '/fs/scratch/PCON0080/yimin/dgscrna'
base <- file.path(root, 'results/hvg_ptc_20260916_v1/ptc_recovery')
stopifnot(file.exists(file.path(base, 'ARCHIVE_STAGED')))
stopifnot(requireNamespace('jsonlite', quietly=TRUE))
dest <- file.path(base, 'inventory_r')
dir.create(dest, recursive=TRUE, showWarnings=FALSE)
packages <- c('Seurat','SeuratObject','Matrix','dbscan','DoubletFinder','future',
              'reticulate','rhdf5','hdf5r','harmony','SeuratDisk','openxlsx',
              'presto','glmGamPoi','jsonlite')
environment <- lapply(packages, function(p) {
  available <- requireNamespace(p, quietly=TRUE)
  list(package=p, available=available,
       version=if (available) as.character(packageVersion(p)) else NA_character_,
       location=if (available) find.package(p) else NA_character_)
})
jsonlite::write_json(list(R=R.version.string, libraries=.libPaths(),
                        job=Sys.getenv('SLURM_JOB_ID'), packages=environment),
                    file.path(dest, 'environment.json'), pretty=TRUE, auto_unbox=TRUE, na='null')
writeLines(capture.output(sessionInfo()), file.path(dest,'sessionInfo.txt'))
if (requireNamespace('DoubletFinder', quietly=TRUE)) {
  writeLines(getNamespaceExports('DoubletFinder'),file.path(dest,'DoubletFinder_exports.txt'))
}

archive <- file.path(base, 'archive')
paths <- list.files(archive, pattern='(full_marker.*\\.RDS|integrated_data\\.RDS)$',
                    recursive=TRUE, full.names=TRUE, ignore.case=TRUE)
inventory <- list()
for (i in seq_along(paths)) {
  path <- paths[[i]]
  key <- sprintf('object_%02d_%s',i,tools::file_path_sans_ext(basename(path)))
  cat(format(Sys.time(),tz='UTC'), 'Reading',path,'\n');flush.console()
  result <- tryCatch({
    obj <- readRDS(path)
    info <- list(source=path,key=key,class=class(obj),length=length(obj),
                 size_bytes=as.numeric(object.size(obj)),names=names(obj))
    writeLines(capture.output(str(obj,max.level=3,list.len=30,give.attr=FALSE)),
               file.path(dest,paste0(key,'.structure.txt')))
    if (inherits(obj,'Seurat')) {
      meta <- methods::slot(obj,'meta.data')
      write.csv(meta,file.path(dest,paste0(key,'.meta.csv')))
      info$n_cells <- ncol(obj)
      info$n_features <- nrow(obj)
      info$active_assay <- methods::slot(obj,'active.assay')
      info$assays <- lapply(methods::slot(obj,'assays'),function(a) {
        list(class=class(a),dim=dim(a),slots=methods::slotNames(a),
             layers=if ('layers' %in% methods::slotNames(a)) names(methods::slot(a,'layers')) else NULL)
      })
      info$reductions <- lapply(methods::slot(obj,'reductions'),function(x) dim(methods::slot(x,'cell.embeddings')))
      info$seurat_version <- as.character(methods::slot(obj,'version'))
      writeLines(capture.output(str(methods::slot(obj,'commands'),max.level=4,list.len=60,give.attr=FALSE)),
                 file.path(dest,paste0(key,'.commands.txt')))
      info$meta_columns <- colnames(meta)
      info$meta_low_cardinality <- lapply(meta[vapply(meta,function(x) length(unique(x))<=200,logical(1))],
                                          function(x) as.list(table(x,useNA='ifany')))
    } else if (is.list(obj)) {
      info$children <- lapply(obj,function(x) list(class=class(x),length=length(x),dim=dim(x),
                         columns=colnames(x),names=if (is.null(dim(x))) names(x) else NULL))
      # Exact marker symbols are exported separately for source comparisons.
      jsonlite::write_json(obj,file.path(dest,paste0(key,'.content.json')),
                          pretty=FALSE,auto_unbox=FALSE,na='null',dataframe='columns')
    }
    rm(obj);gc()
    info
  },error=function(e) list(source=path,key=key,error=conditionMessage(e)))
  inventory[[key]] <- result
  jsonlite::write_json(inventory,file.path(dest,'objects.json'),pretty=TRUE,auto_unbox=TRUE,na='null')
}
writeLines(format(Sys.time(),tz='UTC'),file.path(dest,'COMPLETE'))
cat('Historical R inventory finished; any recorded object errors require review.\n')
