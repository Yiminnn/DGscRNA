# Generalizes research/reference_examples/verify_r_artifacts.R only by matching
# density libraries by name. The original exact object comparison is retained.
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
extra <- Sys.getenv('DGSCRNA_REFERENCE_R_LIB')
if(nzchar(extra)) .libPaths(c(extra,.libPaths()))
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==3L)
left <- readRDS(file.path(args[[1]],'DEG.rds'))
right <- readRDS(file.path(args[[2]],'DEG.rds'))
stopifnot(identical(left,right))
left <- readRDS(file.path(args[[1]],'density_scores.rds'))
right <- readRDS(file.path(args[[2]],'density_scores.rds'))
libraries <- jsonlite::fromJSON(args[[3]])
stopifnot(length(libraries)>0L,!anyDuplicated(libraries),
          setequal(names(left),libraries),all(libraries %in% names(right)))
for(name in libraries) stopifnot(identical(left[[name]],right[[name]]))
cat('R_DEG_AND_NAMED_DENSITY_OBJECTS_EXACT\n')
