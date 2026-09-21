#!/usr/bin/env Rscript
# Rebuild the historical dbscan integer-index fix in the new environment.
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
args <- commandArgs(trailingOnly=TRUE)
tar <- normalizePath(args[[1]], mustWork=TRUE)
work <- args[[2]]
lib <- .libPaths()[[1]]
stopifnot(digest::digest(file=tar,algo='sha256')=='b4eab5a7ec4bdc2b65a1e89ac1cdd5860fe64f407953f19fc8e0af2922f0c9a6')
dir.create(work, recursive=TRUE, showWarnings=FALSE)
untar(tar, exdir=work)
cpp <- file.path(work,'dbscan/src/mst.cpp')
original <- readLines(cpp)
needle <- '      double the_weight = x_dist[LT_POS0(n, node, i)];'
stopifnot(sum(original==needle)==1L)
fixed <- original
fixed[fixed==needle] <- '      double the_weight = x_dist[LT_POS0(n, static_cast<R_xlen_t>(node), static_cast<R_xlen_t>(i))];'
writeLines(fixed,cpp)
description <- file.path(work,'dbscan/DESCRIPTION')
d <- readLines(description);stopifnot(sum(d=='Version: 1.2.6')==1L)
d[d=='Version: 1.2.6'] <- 'Version: 1.2.6.9001';writeLines(d,description)
Sys.setenv(MAKEFLAGS='-j4')
install.packages(file.path(work,'dbscan'),lib=lib,repos=NULL,type='source')
stopifnot(as.character(packageVersion('dbscan'))=='1.2.6.9001')
jsonlite::write_json(list(package='dbscan',version='1.2.6.9001',library=lib,
 official_tar_sha256=digest::digest(file=tar,algo='sha256'),
 patched_cpp_sha256=digest::digest(file=cpp,algo='sha256'),
 change='Both lower-triangle indices cast to R_xlen_t before integer products; original algorithm and parameters unchanged',
 job=Sys.getenv('SLURM_JOB_ID')),file.path(work,'install_manifest.json'),pretty=TRUE,auto_unbox=TRUE)
