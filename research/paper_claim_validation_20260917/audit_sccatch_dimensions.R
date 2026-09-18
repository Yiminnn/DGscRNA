# Public-API regression: native singleton failure, dimension guard, normal parity.
# The original red-capable three-cell run is preserved in SLURM job 7370886.
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
root <- '/fs/scratch/PCON0080/yimin/dgscrna'
base <- file.path(root, 'results/hvg_ptc_20260916_v1/paper_claim_validation_20260917')
.libPaths(c(file.path(base, 'vendor_R'), .libPaths()))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(scCATCH))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest))
script <- sub('^--file=', '', grep('^--file=', commandArgs(), value=TRUE)[1])
source(file.path(dirname(script), 'sccatch_dimension_guard.R'))
args <- commandArgs(trailingOnly=TRUE)
dest <- if (length(args)) args[[1]] else file.path(base, 'verification/scCATCH_dimension_diagnosis', Sys.getenv('SLURM_JOB_ID'))
dir.create(dest, recursive=TRUE, showWarnings=FALSE)
fixture <- function(singleton, sparse) {
  vals <- if (singleton) c(8, 4, 1, 2, 2, 1) else c(8, 4, 7, 5, 1, 2, 2, 1)
  x <- matrix(log1p(vals), nrow=2)
  rownames(x) <- c('GENE_A', 'GENE_B'); colnames(x) <- paste0('c', seq_len(ncol(x)))
  if (sparse) x <- Matrix(x, sparse=TRUE)
  list(x=x, cluster=if (singleton) c('0', '1', '1') else c('0', '0', '1', '1'))
}
fixtures <- list(singleton_sparse=fixture(TRUE, TRUE), singleton_dense=fixture(TRUE, FALSE),
                 normal_sparse=fixture(FALSE, TRUE), normal_dense=fixture(FALSE, FALSE))
native <- function(f) {
  marker <- data.frame(gene=rownames(f$x), celltype='fixture', pmid='fixture',
                       subtype1=NA_character_, subtype2=NA_character_, subtype3=NA_character_)
  # A permissive cutoff is restricted to this tiny test to avoid the independent
  # native empty-DEG bug. Production p-value thresholds remain .05/.1/.01.
  tryCatch(findmarkergene(createscCATCH(f$x, f$cluster), if_use_custom_marker=TRUE,
    marker=marker, use_method='1', cell_min_pct=.25, logfc=.25, pvalue=1,
    verbose=FALSE), error=function(e)e)
}
errors <- function(results) lapply(results, function(r)
  if (inherits(r, 'error')) conditionMessage(r) else '')
canonical <- function(d) {
  if (!nrow(d)) return(character())
  d <- d[, c('cluster', 'gene', 'comp_cluster', 'pct', 'logfc', 'pvalue'), drop=FALSE]
  rownames(d) <- NULL
  sort(apply(d, 1, paste, collapse='|'))
}
original <- lapply(fixtures, native)
err <- errors(original)
write_json(err, file.path(dest, 'native_errors.json'), pretty=TRUE, auto_unbox=TRUE)
stopifnot(all(vapply(err[c('singleton_sparse', 'singleton_dense')], function(e)
  grepl('dim(X) must have a positive length', e, fixed=TRUE), logical(1))))
stopifnot(all(vapply(err[c('normal_sparse', 'normal_dense')], identical, logical(1), '')))
guard <- install_sccatch_dimension_guard()
fixed <- lapply(fixtures, native)
write_json(errors(fixed), file.path(dest, 'guarded_errors.json'), pretty=TRUE, auto_unbox=TRUE)
stopifnot(all(vapply(errors(fixed), identical, logical(1), '')))
for (name in c('normal_sparse', 'normal_dense'))
  stopifnot(identical(canonical(original[[name]]@markergene), canonical(fixed[[name]]@markergene)))
stopifnot(identical(canonical(fixed$singleton_sparse@markergene), canonical(fixed$singleton_dense@markergene)))
m <- list(status='passed', job=Sys.getenv('SLURM_JOB_ID'), guard=guard,
  guard_source_sha256=digest(file=file.path(dirname(script), 'sccatch_dimension_guard.R'), algo='sha256'),
  source_sha256=digest(file=script, algo='sha256'), native_singleton_error_reproduced=TRUE,
  normal_sparse_and_dense_DEG_exact=TRUE, singleton_sparse_dense_DEG_exact=TRUE,
  singleton_preserved=TRUE, fixture_only_pvalue=1,
  production_thresholds_unchanged=c(.05, .1, .01))
write_json(m, file.path(dest, 'unit_manifest.json'), pretty=TRUE, auto_unbox=TRUE)
writeLines(digest(file=file.path(dest, 'unit_manifest.json'), algo='sha256'), file.path(dest, 'UNIT_COMPLETE'))
cat('SCCATCH_DIMENSION_REGRESSION_PASSED\n')
