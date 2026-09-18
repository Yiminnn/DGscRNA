# Official SCINA 1.2.0 on the same full normalized RNA and frozen panel genes.
# Empty/constant or overlapping signatures are recorded, never replaced by truth.
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest))
root<-'/fs/scratch/PCON0080/yimin/dgscrna'
base<-file.path(root,'results/hvg_ptc_20260916_v1/paper_claim_validation_20260917')
.libPaths(c(file.path(base,'vendor_R'),.libPaths()))
suppressPackageStartupMessages(library(SCINA))
script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
source(file.path(dirname(script),'scina_numerics.R'))
guard_manifest<-fromJSON(file.path(base,'verification/SCINA_stability_v3/manifest.json'),simplifyVector=FALSE)
stopifnot(guard_manifest$status=='passed',
  guard_manifest$numerical_source_sha256==digest(file=file.path(dirname(script),'scina_numerics.R'),algo='sha256'))
guarded_SCINA<-make_numerically_guarded_SCINA(base)
a<-commandArgs(trailingOnly=TRUE);sample<-a[[1]]
index<-if(length(a)>1L)as.integer(a[[2]]) else as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
libs<-fromJSON(file.path(base,'markers/libraries.json'),simplifyVector=FALSE)
lib<-names(libs)[index+1L];signatures<-lapply(libs[[index+1L]],unlist,use.names=FALSE)
dest<-file.path(base,'comparators/SCINA',sample,sprintf('L%02d',index))
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
if(file.exists(file.path(dest,'COMPLETE')))quit(status=0)
previous<-list.files(dest,full.names=TRUE);previous<-previous[!dir.exists(previous)]
if(length(previous)) {
  archive<-file.path(dest,'attempts',paste0('before_',Sys.getenv('SLURM_JOB_ID')))
  dir.create(archive,recursive=TRUE,showWarnings=FALSE);stopifnot(all(file.copy(previous,archive)))
}
obj<-readRDS(file.path(base,'GBM',sample,'hvg2000/expression_PCA30.rds'))
genes<-intersect(rownames(obj),unique(unlist(signatures,use.names=FALSE)))
X<-as.matrix(obj[['RNA']]@data[genes,,drop=FALSE]);rm(obj);gc()
pred<-data.frame(cell_id=colnames(X),stringsAsFactors=FALSE);arms<-list()
for(overlap in c(1L,0L)) {
  key<-paste0('overlap',overlap);started<-Sys.time()
  logfile<-file.path(dest,paste0(key,'.log'));warnings<-character()
  # Run the package's own preprocessing; remove its resulting zero-length
  # signatures explicitly because 1.2.0 otherwise reaches chol on a 0x0 matrix.
  q<-SCINA:::check.inputs(X,signatures,100L,10L,.99,1,overlap,logfile)
  stopifnot(q$qual==1)
  sig<-q$sig[lengths(q$sig)>0L]
  status<-'completed';error<-NULL;native_error<-NULL;solver<-'official_unmodified';labels<-rep('Unknown',ncol(X))
  probabilities<-NULL
  if(!length(sig))status<-'no_supported_nonconstant_signatures' else {
    fit<-tryCatch(withCallingHandlers(SCINA::SCINA(X,sig,max_iter=100L,convergence_n=10L,
      convergence_rate=.99,sensitivity_cutoff=1,rm_overlap=overlap,allow_unknown=1,
      log_file=logfile),warning=function(w){warnings<<-c(warnings,conditionMessage(w))}),
      error=function(e)e)
    if(inherits(fit,'error')) {
      native_error<-conditionMessage(fit);solver<-'audited_numerical_boundary_guard'
      stopifnot(file.copy(logfile,paste0(logfile,'.native_failure'),overwrite=TRUE))
      fit<-tryCatch(guarded_SCINA(X,sig,max_iter=100L,convergence_n=10L,convergence_rate=.99,
        sensitivity_cutoff=1,rm_overlap=overlap,allow_unknown=1,log_file=logfile),error=function(e)e)
    }
    if(inherits(fit,'error')){status<-'implementation_error';error<-conditionMessage(fit)} else {
      stopifnot(length(fit$cell_labels)==ncol(X),all(is.finite(fit$probabilities)))
      labels<-fit$cell_labels;labels[labels=='unknown']<-'Unknown'
      probabilities<-fit$probabilities
      saveRDS(list(probabilities=probabilities,cell_id=colnames(X),signatures=sig),
        file.path(dest,paste0(key,'_probabilities.rds')),compress=FALSE)
    }
  }
  pred[[key]]<-labels
  arms[[key]]<-list(method='SCINA',library=lib,budget='full_RNA',route='cell_level',
    cutoff=key,status=status,error=error,native_error=native_error,solver=solver,n_input_panels=length(signatures),
    n_supported_panels=length(sig),dropped_panels=setdiff(names(signatures),names(sig)),
    rm_overlap=overlap,allow_unknown=1,max_iter=100L,convergence_n=10L,convergence_rate=.99,
    sensitivity_cutoff=1,warnings=unique(warnings),elapsed_seconds=as.numeric(difftime(Sys.time(),started,units='secs')))
  write_json(arms,file.path(dest,'arm_progress.json'),pretty=TRUE,auto_unbox=TRUE)
  cat('SCINA_ARM',sample,lib,key,status,'\n');flush.console()
}
con<-gzfile(file.path(dest,'predictions.csv.gz'),'wt');write.csv(pred,con,row.names=FALSE);close(con)
m<-list(status=if(any(vapply(arms,function(x)x$status=='implementation_error',logical(1))))'needs_review' else 'completed',
  method='SCINA',version=as.character(packageVersion('SCINA')),sample=sample,library=lib,arms=arms,
  n_cells=ncol(X),input='Seurat RNA log1p library-size normalized to 10000; shared frozen signature universe',
  preprocessing='Official check.inputs plus dropping resulting empty signatures; no truth-based gene filtering',
  numerical_guard_audit_sha256=digest(file=file.path(base,'verification/SCINA_stability_v3/manifest.json'),algo='sha256'),
  marker_sha256=digest(file=file.path(base,'markers/libraries.json'),algo='sha256'),
  predictions_sha256=digest(file=file.path(dest,'predictions.csv.gz'),algo='sha256'),
  job=Sys.getenv('SLURM_JOB_ID'),source_sha256=digest(file=sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]),algo='sha256'))
write_json(m,file.path(dest,'manifest.json'),pretty=TRUE,auto_unbox=TRUE)
if(m$status=='completed')writeLines(digest(file=file.path(dest,'manifest.json'),algo='sha256'),file.path(dest,'COMPLETE'))
if(m$status!='completed')stop('SCINA implementation error requires review; not a valid zero-performance result')
