stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
base <- '/fs/scratch/PCON0080/yimin/dgscrna/results/hvg_ptc_20260916_v1/ptc_experiments'
mode <- commandArgs(trailingOnly=TRUE)[[1]]
lib <- file.path(base,if(mode=='official')'vendor_R' else 'vendor_R_dbscan64')
.libPaths(c(lib,.libPaths()))
library(dbscan)
stopifnot(as.character(packageVersion('dbscan'))==if(mode=='official')'1.2.6' else '1.2.6.9001')
dest <- file.path(base,'verification/dbscan_index64');dir.create(dest,recursive=TRUE,showWarnings=FALSE)
set.seed(712)
x <- rbind(matrix(rnorm(1200,0,.15),ncol=3),matrix(rnorm(1200,4,.3),ncol=3),matrix(runif(150,-1,5),ncol=3))
fit <- hdbscan(x,minPts=20)
if(mode=='official') {
  saveRDS(fit,file.path(dest,'official_small_reference.rds'))
} else {
  reference <- readRDS(file.path(dest,'official_small_reference.rds'))
  for(n in c('cluster','coredist','cluster_scores','membership_prob','outlier_scores'))
    stopifnot(isTRUE(all.equal(fit[[n]],reference[[n]],tolerance=0)))
  # Compile the exact patched index expression and check against double-precision
  # integer arithmetic across the overflow boundary, without allocating an N^2 matrix.
  Rcpp::cppFunction('Rcpp::NumericVector index64(Rcpp::NumericVector ns,Rcpp::NumericVector is,Rcpp::NumericVector js) {
    Rcpp::NumericVector result(ns.size());
    for (int k=0;k<ns.size();++k) {
      R_xlen_t n=ns[k], i=is[k], j=js[k];
      result[k]=(i)==(j) ? 0 : (i)<(j) ? (n)*(i)-((i)+1)*(i)/2+(j)-(i)-1 : (n)*(j)-((j)+1)*(j)/2+(i)-(j)-1;
    }
    return result;
  }')
  n <- rep(48255,100004);i <- c(0,46340,46341,48253,sample(0:48254,100000,replace=TRUE))
  j <- c(48254,48254,48254,48254,sample(0:48254,100000,replace=TRUE))
  lo <- pmin(as.double(i),as.double(j));hi <- pmax(as.double(i),as.double(j))
  expected <- ifelse(i==j,0,n*lo-(lo+1)*lo/2+hi-lo-1)
  stopifnot(identical(index64(n,i,j),expected),all(expected>=0),all(expected<n*(n-1)/2))
  jsonlite::write_json(list(status='passed',official_version='1.2.6',patched_version='1.2.6.9001',
    small_clustering_fields_exact=c('cluster','coredist','cluster_scores','membership_prob','outlier_scores'),
    n_index_checks=length(n),overflow_boundary_covered=TRUE,
    original_integer_product_at46341=as.double(46341)*46342,
    int32_max=.Machine$integer.max,job=Sys.getenv('SLURM_JOB_ID')),
    file.path(dest,'verification.json'),auto_unbox=TRUE,pretty=TRUE)
  writeLines(digest::digest(file=file.path(dest,'verification.json'),algo='sha256'),file.path(dest,'COMPLETE'))
}
