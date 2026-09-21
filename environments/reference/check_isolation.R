args <- commandArgs(trailingOnly=TRUE)
prefix <- normalizePath(args[[1]],mustWork=TRUE)
packages <- c('Seurat','SeuratObject','Matrix','future','jsonlite','digest','uwot','limma','dbscan')
for(p in packages)suppressPackageStartupMessages(library(p,character.only=TRUE))
namespaces <- loadedNamespaces()
locations <- vapply(namespaces,function(p)normalizePath(find.package(p)),character(1))
stopifnot(all(startsWith(locations,paste0(prefix,'/'))))
versions <- vapply(namespaces,function(p)as.character(packageVersion(p)),character(1))
jsonlite::write_json(list(status='all_loaded_R_packages_in_new_prefix',prefix=prefix,
 libraries=.libPaths(),loaded_namespace_paths=as.list(locations),versions=as.list(versions),
 R_version=R.version.string,job=Sys.getenv('SLURM_JOB_ID')),args[[2]],auto_unbox=TRUE,pretty=TRUE)
