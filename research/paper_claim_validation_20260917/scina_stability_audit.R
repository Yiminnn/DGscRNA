stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
suppressPackageStartupMessages(library(Seurat));suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest))
root<-'/fs/scratch/PCON0080/yimin/dgscrna';base<-file.path(root,'results/hvg_ptc_20260916_v1/paper_claim_validation_20260917')
.libPaths(c(file.path(base,'vendor_R'),.libPaths()));suppressPackageStartupMessages(library(SCINA))
dest<-file.path(base,'verification/SCINA_stability_v3');dir.create(dest,recursive=TRUE,showWarnings=FALSE)
obj<-readRDS(file.path(base,'GBM/TKU4163/hvg2000/expression_PCA30.rds'))
libs<-fromJSON(file.path(base,'markers/libraries.json'),simplifyVector=FALSE)
stable_density_ratio<-function(e,mu1,mu2,inverse_sigma1,inverse_sigma2) {
  original<-SCINA:::density_ratio(e,mu1,mu2,inverse_sigma1,inverse_sigma2)
  if(all(is.finite(original)))return(original)
  a<-colSums((e-mu1)*(inverse_sigma1%*%(e-mu1)))
  b<-colSums((e-mu2)*(inverse_sigma2%*%(e-mu2)))
  # The original solver ties the two covariance matrices at every update.
  # Their equal log-determinants cancel exactly; evaluating det() separately
  # can overflow to Inf and produce Inf-Inf for large marker signatures.
  if(identical(inverse_sigma1,inverse_sigma2))log_det_ratio<-0 else
    log_det_ratio<-as.numeric(determinant(inverse_sigma2,logarithm=TRUE)$modulus-
                              determinant(inverse_sigma1,logarithm=TRUE)$modulus)
  logratio<--.5*(a-b+log_det_ratio)
  fixed<-exp(pmax(log(1e-200),pmin(log(1e200),logratio)))
  original[!is.finite(original)]<-fixed[!is.finite(original)]
  original
}
script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
source(file.path(dirname(script),'scina_numerics.R'))
fn<-make_numerically_guarded_SCINA(base)
checks<-list()
for(i in 0:15)for(overlap in c(1L,0L)) {
  signatures<-lapply(libs[[i+1L]],unlist,use.names=FALSE)
  X<-as.matrix(obj[['RNA']]@data[intersect(rownames(obj),unique(unlist(signatures))),,drop=FALSE])
  key<-sprintf('L%02d_overlap%d',i,overlap);log<-file.path(dest,paste0(key,'.log'))
  q<-SCINA:::check.inputs(X,signatures,100L,10L,.99,1,overlap,log);sig<-q$sig[lengths(q$sig)>0L]
  failure<-NULL
  fit<-tryCatch(withCallingHandlers(fn(X,sig,rm_overlap=overlap,log_file=log),error=function(e) {
    for(f in sys.frames())if(exists('theta',envir=f,inherits=FALSE)&&exists('prob_mat',envir=f,inherits=FALSE)) {
      pm<-get('prob_mat',envir=f);th<-get('theta',envir=f)
      failure<<-list(message=conditionMessage(e),iteration=get('iter',envir=f),
        probability_nonfinite=sum(!is.finite(pm)),zero_probability_rows=which(rowSums(pm)==0),
        nonfinite_theta=which(vapply(th,function(v)any(!is.finite(v$mean)),logical(1))))
    }
  }),error=function(e)e)
  original<-file.path(base,'comparators/SCINA/TKU4163',sprintf('L%02d',i),paste0('overlap',overlap,'_probabilities.rds'))
  if(inherits(fit,'error'))check<-list(status='error',message=conditionMessage(fit),details=failure) else {
    check<-list(status='completed',finite=all(is.finite(fit$probabilities)))
    if(file.exists(original)) {
      ref<-readRDS(original);stopifnot(identical(dimnames(ref$probabilities),dimnames(fit$probabilities)))
      check$max_abs_delta<-max(abs(ref$probabilities-fit$probabilities))
      check$probability_parity<-check$max_abs_delta<1e-8
      # Record all comparisons before applying the final release gate.
    }
    saveRDS(fit,file.path(dest,paste0(key,'.rds')))
  }
  checks[[key]]<-check;write_json(checks,file.path(dest,'checks.json'),pretty=TRUE,auto_unbox=TRUE)
  cat('STABILITY_AUDIT',key,check$status,'\n');flush.console()
}
write_json(list(status=if(all(vapply(checks,function(v)v$status=='completed' && (is.null(v$probability_parity)||v$probability_parity),logical(1))))'passed' else 'needs_review',
  checks=checks,job=Sys.getenv('SLURM_JOB_ID'),source_sha256=digest(file=script,algo='sha256'),
  numerical_source_sha256=digest(file=file.path(dirname(script),'scina_numerics.R'),algo='sha256')),
  file.path(dest,'manifest.json'),pretty=TRUE,auto_unbox=TRUE)
