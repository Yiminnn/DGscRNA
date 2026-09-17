stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
base <- '/fs/scratch/PCON0080/yimin/dgscrna/results/hvg_ptc_20260916_v1/ptc_experiments'
tar <- file.path(base,'source_snapshot/dbscan_1.2.6.tar.gz')
stopifnot(digest::digest(file=tar,algo='sha256')=='b4eab5a7ec4bdc2b65a1e89ac1cdd5860fe64f407953f19fc8e0af2922f0c9a6')
src <- file.path(base,'source_snapshot/dbscan_index64_source');dir.create(src,showWarnings=FALSE)
untar(tar,exdir=src)
cpp <- file.path(src,'dbscan/src/mst.cpp')
original <- readLines(cpp)
needle <- '      double the_weight = x_dist[LT_POS0(n, node, i)];'
stopifnot(sum(original==needle)==1L)
fixed <- original
fixed[fixed==needle] <- '      double the_weight = x_dist[LT_POS0(n, static_cast<R_xlen_t>(node), static_cast<R_xlen_t>(i))];'
writeLines(fixed,cpp)
description <- file.path(src,'dbscan/DESCRIPTION')
d <- readLines(description);stopifnot(sum(d=='Version: 1.2.6')==1L)
d[d=='Version: 1.2.6'] <- 'Version: 1.2.6.9001';writeLines(d,description)
lib <- file.path(base,'vendor_R_dbscan64');dir.create(lib,showWarnings=FALSE)
Sys.setenv(MAKEFLAGS='-j8')
install.packages(file.path(src,'dbscan'),lib=lib,repos=NULL,type='source')
.libPaths(c(lib,.libPaths()))
stopifnot(as.character(packageVersion('dbscan'))=='1.2.6.9001')
jsonlite::write_json(list(package='dbscan',version='1.2.6.9001',library=lib,
  official_tar_sha256=digest::digest(file=tar,algo='sha256'),
  patched_cpp_sha256=digest::digest(file=cpp,algo='sha256'),
  change='Cast both MST lower-triangle indices to R_xlen_t before LT_POS0 products; no algorithm/parameter change',
  mechanism='MST node and loop index were int; ((i)+1)*(i) overflows32-bit before promotion even when n is R_xlen_t',
  affected_boundary='i>=46341; MTN48255 cells exceeds the boundary, TUT44149 does not',
  diff=list(before=needle,after=fixed[original==needle]),
  source='https://github.com/mhahsler/dbscan/blob/dbscan_1.2.6/src/mst.cpp',
  header='https://github.com/mhahsler/dbscan/blob/dbscan_1.2.6/src/lt.h',
  job=Sys.getenv('SLURM_JOB_ID')),file.path(base,'source_snapshot/dbscan_index64_install.json'),
  pretty=TRUE,auto_unbox=TRUE)
writeLines('1.2.6.9001',file.path(lib,'DBSCAN_INSTALLED'))
