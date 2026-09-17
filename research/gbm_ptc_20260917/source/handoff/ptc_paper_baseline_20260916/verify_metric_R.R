stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
dest <- '/fs/scratch/PCON0080/yimin/dgscrna/results/hvg_ptc_20260916_v1/ptc_paper_baseline'
source('/fs/scratch/PCON0080/yimin/dgscrna/handoff/ptc_paper_baseline_20260916/reference_sources/MLmetrics_1.1.1_Classification.R')
writeLines(c('MLmetrics 1.1.1 archived CRAN source; no package installation or edited function body',deparse(F1_Score),
             deparse(Precision),deparse(Recall)),file.path(dest,'MLmetrics_actual_definitions.txt'))
d <- read.csv(gzfile(file.path(dest,'original_DG_binary_pairs_for_R.csv.gz')))
scopes <- c('Overall',unique(d$scope))
res <- do.call(rbind,lapply(scopes,function(scope) {
  x <- if(scope=='Overall')d else d[d$scope==scope,]
  data.frame(scope=scope,MLmetrics_default=F1_Score(x$truth,x$prediction),
    MLmetrics_positive0=F1_Score(x$truth,x$prediction,positive='0'),
    MLmetrics_positive1=F1_Score(x$truth,x$prediction,positive='1'))
}))
write.csv(res,file.path(dest,'MLmetrics_original_DG_verification.csv'),row.names=FALSE)
print(res)
