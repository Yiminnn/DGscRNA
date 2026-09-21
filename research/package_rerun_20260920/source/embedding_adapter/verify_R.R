#!/usr/bin/env Rscript
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
.libPaths(.Library,include.site=FALSE)
args<-commandArgs(trailingOnly=TRUE)
pairs<-jsonlite::fromJSON(args[[1]],simplifyVector=FALSE)
rows<-list()
for(pair in pairs) {
  for(name in c('DEG.rds','density_scores.rds')) {
    left<-readRDS(file.path(pair$actual,name))
    right<-readRDS(file.path(pair$expected,name))
    stopifnot(identical(left,right))
  }
  rows[[length(rows)+1L]]<-list(space=pair$space,partition=pair$partition,
      DEG_exact=TRUE,density_exact=TRUE,
      actual_DEG_sha256=digest::digest(file=file.path(pair$actual,'DEG.rds'),algo='sha256'),
      actual_density_sha256=digest::digest(file=file.path(pair$actual,'density_scores.rds'),algo='sha256'),
      expected_DEG_sha256=digest::digest(file=file.path(pair$expected,'DEG.rds'),algo='sha256'),
      expected_density_sha256=digest::digest(file=file.path(pair$expected,'density_scores.rds'),algo='sha256'))
}
jsonlite::write_json(rows,args[[2]],pretty=TRUE,auto_unbox=TRUE)
cat('R_DEG_DENSITY_EXACT',length(rows),'partitions\n')
