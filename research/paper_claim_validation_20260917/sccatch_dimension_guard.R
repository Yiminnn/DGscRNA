# Process-local compatibility fix for the confirmed scCATCH 3.2.2 singleton bug.
# The two column selections retain matrix dimensions. No algorithm is replaced.
install_sccatch_dimension_guard <- function(expected_original=NULL) {
  stopifnot(as.character(packageVersion('scCATCH')) == '3.2.2')
  ns <- asNamespace('scCATCH')
  edits <- list(
    .get_marker_scCATCH1=c(
      'ndata1 <- ndata[, meta[meta$cluster == clu_num[i], ]$cell]',
      'ndata1 <- ndata[, meta[meta$cluster == clu_num[i], ]$cell, drop = FALSE]'),
    .get_marker=c(
      'ndata2 <- ndata_temp[, meta[meta$cluster == clu_pair1$cluster2[j], ]$cell]',
      'ndata2 <- ndata_temp[, meta[meta$cluster == clu_pair1$cluster2[j], ]$cell, drop = FALSE]'))
  records <- list()
  for (name in names(edits)) {
    fun <- get(name, envir=ns)
    original <- paste(deparse(body(fun), width.cutoff=500L), collapse='\n')
    before <- digest::digest(original, algo='sha256', serialize=FALSE)
    if (!is.null(expected_original)) stopifnot(identical(before, expected_original[[name]]))
    edit <- edits[[name]]
    stopifnot(sum(gregexpr(edit[[1]], original, fixed=TRUE)[[1]] > 0L) == 1L)
    patched <- sub(edit[[1]], edit[[2]], original, fixed=TRUE)
    body(fun) <- parse(text=patched)[[1]]
    unlockBinding(name, ns)
    assign(name, fun, envir=ns)
    lockBinding(name, ns)
    records[[name]] <- list(original_sha256=before,
      patched_sha256=digest::digest(paste(deparse(body(fun), width.cutoff=500L), collapse='\n'),
                                   algo='sha256', serialize=FALSE),
      original_expression=edit[[1]], guarded_expression=edit[[2]])
  }
  list(kind='scCATCH_3.2.2_singleton_dimension_preservation', functions=records,
    retained='All cells, original clusters, markers, thresholds, tests and native scoring',
    scope='Runtime namespace of this R process only; installed package files unchanged')
}
