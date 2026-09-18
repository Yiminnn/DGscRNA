# Exact parity with the completed native TKU4163 public-function audit.
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest))
a <- commandArgs(trailingOnly=TRUE)
original <- a[[1]]; guarded <- a[[2]]
om <- fromJSON(file.path(original, 'audit_manifest.json'), simplifyVector=FALSE)
gm <- fromJSON(file.path(guarded, 'cohort_manifest.json'), simplifyVector=FALSE)
stopifnot(identical(om$sample, gm$sample), identical(om$n_cells, gm$n_cells),
          length(om$arms) == 48L, length(gm$arms) == 48L,
          identical(om$arms, gm$arms), identical(om$thresholds, gm$thresholds))
op <- file.path(original, 'audit_predictions.csv.gz')
gp <- file.path(guarded, 'cohort_predictions.csv.gz')
stopifnot(identical(digest(file=op, algo='sha256'), om$predictions_sha256),
          identical(digest(file=gp, algo='sha256'), gm$predictions_sha256))
old <- read.csv(gzfile(op), colClasses='character', check.names=FALSE)
new <- read.csv(gzfile(gp), colClasses='character', check.names=FALSE)
stopifnot(identical(old, new), ncol(old) == 49L, nrow(old) == 178L)
canonical <- function(d) {
  if (!nrow(d)) return(character())
  d <- d[, c('cluster', 'gene', 'comp_cluster', 'pct', 'logfc', 'pvalue'), drop=FALSE]
  rownames(d) <- NULL
  sort(apply(d, 1, paste, collapse='|'))
}
od <- file.path(original, 'official_union_DEG.rds')
gd <- file.path(guarded, 'official_union_DEG.rds')
stopifnot(identical(canonical(readRDS(od)), canonical(readRDS(gd))))
files <- c(op, gp, od, gd, file.path(original, 'audit_manifest.json'),
           file.path(guarded, 'cohort_manifest.json'))
write_json(list(status='passed', job=Sys.getenv('SLURM_JOB_ID'), sample=om$sample,
  n_cells=nrow(old), n_conditions=48L, exact_cell_order=TRUE,
  exact_all_native_labels=TRUE, exact_all_arms=TRUE, exact_union_DEG=TRUE,
  files=setNames(lapply(files, function(p)digest(file=p, algo='sha256')), files)),
  file.path(guarded, 'real_parity_manifest.json'), pretty=TRUE, auto_unbox=TRUE)
cat('SCCATCH_REAL_48_ARM_EXACT_PARITY_PASSED\n')
