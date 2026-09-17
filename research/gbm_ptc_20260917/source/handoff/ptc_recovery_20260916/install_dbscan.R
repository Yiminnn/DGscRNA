stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
base <- '/fs/scratch/PCON0080/yimin/dgscrna/results/hvg_ptc_20260916_v1/ptc_experiments'
lib <- file.path(base,'vendor_R');dir.create(lib,recursive=TRUE,showWarnings=FALSE)
src <- file.path(base,'source_snapshot/dbscan_1.2.6.tar.gz')
url <- 'https://cran.r-project.org/src/contrib/dbscan_1.2.6.tar.gz'
download.file(url,src,mode='wb',method='libcurl')
Sys.setenv(MAKEFLAGS='-j8')
install.packages(src,lib=lib,repos=NULL,type='source')
.libPaths(c(lib,.libPaths()))
stopifnot(as.character(packageVersion('dbscan'))=='1.2.6')
jsonlite::write_json(list(package='dbscan',version='1.2.6',library=lib,url=url,
  tar_sha256=digest::digest(file=src,algo='sha256'),
  reason='Upstream fix for MST root-node out-of-bounds write in1.2.5',
  upstream='https://github.com/mhahsler/dbscan/releases/tag/dbscan_1.2.6',
  job=Sys.getenv('SLURM_JOB_ID')),file.path(base,'source_snapshot/dbscan_install.json'),
  pretty=TRUE,auto_unbox=TRUE)
writeLines('1.2.6',file.path(lib,'DBSCAN_INSTALLED'))
