stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
args <- commandArgs(trailingOnly=TRUE)
for(name in c('DEG.rds','density_scores.rds')) {
  left <- readRDS(file.path(args[[1]],name))
  right <- readRDS(file.path(args[[2]],name))
  stopifnot(identical(left,right))
}
cat('R_DEG_AND_DENSITY_OBJECTS_EXACT\n')
