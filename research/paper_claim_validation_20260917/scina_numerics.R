# SCINA1.2.0 numerical boundary repair. No package/global namespace is modified.
# 1. Preserve original density ratios whenever finite; only nonfinite ratios
#    caused by equal covariance determinant overflow use exact cancellation.
# 2. At exactly zero posterior mass, the component mean is unidentifiable;
#    retain its preceding value instead of calculating 0/0. Positive-mass
#    updates, covariance updates, probabilities, stopping rules stay original.
make_numerically_guarded_SCINA<-function(base) {
  src<-file.path(base,'vendor_build/SCINA_1.2.0/SCINA/R/EM_model.R')
  txt<-paste(readLines(src),collapse='\n')
  old1<-'theta[[i]]$mean[,1]=(exp[signatures[[i]],]%*%prob_mat[i,])/sum(prob_mat[i,])'
  old2<-'theta[[i]]$mean[,2]=(exp[signatures[[i]],]%*%(1-prob_mat[i,]))/sum(1-prob_mat[i,])'
  stopifnot(length(gregexpr(old1,txt,fixed=TRUE)[[1]])==1L,grepl(old1,txt,fixed=TRUE),
            length(gregexpr(old2,txt,fixed=TRUE)[[1]])==1L,grepl(old2,txt,fixed=TRUE))
  txt<-sub(old1,paste0('if(sum(prob_mat[i,])>0) ',old1),txt,fixed=TRUE)
  txt<-sub(old2,paste0('if(sum(1-prob_mat[i,])>0) ',old2),txt,fixed=TRUE)
  env<-new.env(parent=asNamespace('SCINA'));eval(parse(text=txt),envir=env)
  env$density_ratio<-function(e,mu1,mu2,inverse_sigma1,inverse_sigma2) {
    original<-SCINA:::density_ratio(e,mu1,mu2,inverse_sigma1,inverse_sigma2)
    if(all(is.finite(original)))return(original)
    a<-colSums((e-mu1)*(inverse_sigma1%*%(e-mu1)))
    b<-colSums((e-mu2)*(inverse_sigma2%*%(e-mu2)))
    if(identical(inverse_sigma1,inverse_sigma2))log_det_ratio<-0 else
      log_det_ratio<-as.numeric(determinant(inverse_sigma2,logarithm=TRUE)$modulus-
                                determinant(inverse_sigma1,logarithm=TRUE)$modulus)
    logratio<--.5*(a-b+log_det_ratio)
    fixed<-exp(pmax(log(1e-200),pmin(log(1e200),logratio)))
    original[!is.finite(original)]<-fixed[!is.finite(original)]
    original
  }
  env$SCINA
}
