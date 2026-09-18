# Cold-input SCINA resource measurement with the frozen glioma marker library.
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
root<-'/fs/scratch/PCON0080/yimin/dgscrna'
base<-file.path(root,'results/hvg_ptc_20260916_v1/paper_claim_validation_20260917')
.libPaths(c(file.path(base,'vendor_R'),.libPaths()))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest))
suppressPackageStartupMessages(library(SCINA))
options(Seurat.object.assay.version='v3')
script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
source(file.path(dirname(script),'scina_numerics.R'))
proof<-fromJSON(file.path(base,'verification/SCINA_stability_v3/manifest.json'))
stopifnot(proof$status=='passed',proof$numerical_source_sha256==digest(file=file.path(dirname(script),'scina_numerics.R'),algo='sha256'))
guarded<-make_numerically_guarded_SCINA(base)
a<-commandArgs(trailingOnly=TRUE);n<-as.integer(a[[1]])
sample<-if(length(a)>1L)a[[2]] else paste0('SCALE_',n)
dest<-if(length(a)>2L)a[[3]] else file.path(base,'scalability',n,'SCINA')
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
inp<-file.path(base,'inputs',sample);im<-fromJSON(file.path(inp,'input_manifest.json'))
stopifnot(im$n_cells==n,readLines(file.path(inp,'INPUT_COMPLETE'))==digest(file=file.path(inp,'input_manifest.json'),algo='sha256'))
cells<-read.csv(file.path(inp,'cells_fit.csv'),stringsAsFactors=FALSE)$cell_id
genes<-read.csv(file.path(inp,'genes.csv'),stringsAsFactors=FALSE)$gene
read_bin<-function(name,what,n,size) {con<-file(file.path(inp,name),'rb');on.exit(close(con));readBin(con,what,n=n,size=size,endian='little')}
counts<-new('dgCMatrix',x=read_bin('x.bin','numeric',im$nnz,8),i=read_bin('i.bin','integer',im$nnz,4),
  p=read_bin('p.bin','integer',n+1L,4),Dim=as.integer(c(im$n_genes,n)),Dimnames=list(genes,cells))
obj<-CreateSeuratObject(counts,min.cells=0,min.features=0)
obj<-NormalizeData(obj,normalization.method='LogNormalize',scale.factor=10000,verbose=FALSE)
libs<-fromJSON(file.path(base,'markers/libraries.json'),simplifyVector=FALSE)
sig<-lapply(libs$CM2_glioma_other,unlist,use.names=FALSE)
genes<-intersect(rownames(obj),unique(unlist(sig,use.names=FALSE)))
X<-as.matrix(obj[['RNA']]@data[genes,,drop=FALSE]);rm(obj,counts);gc()
logfile<-file.path(dest,'official_SCINA.log')
q<-SCINA:::check.inputs(X,sig,100L,10L,.99,1,1L,logfile);stopifnot(q$qual==1)
sig<-q$sig[lengths(q$sig)>0L];stopifnot(length(sig)>0L)
solver<-'official_unmodified';native_error<-NULL
fit<-tryCatch(SCINA::SCINA(X,sig,max_iter=100L,convergence_n=10L,convergence_rate=.99,
    sensitivity_cutoff=1,rm_overlap=1L,allow_unknown=1L,log_file=logfile),error=function(e)e)
if(inherits(fit,'error')) {
  native_error<-conditionMessage(fit);solver<-'audited_numerical_boundary_guard'
  stopifnot(file.copy(logfile,paste0(logfile,'.native_failure'),overwrite=TRUE))
  fit<-guarded(X,sig,max_iter=100L,convergence_n=10L,convergence_rate=.99,
    sensitivity_cutoff=1,rm_overlap=1L,allow_unknown=1L,log_file=logfile)
}
stopifnot(length(fit$cell_labels)==n,all(is.finite(fit$probabilities)))
con<-gzfile(file.path(dest,'predictions.csv.gz'),'wt')
write.csv(data.frame(cell_id=cells,prediction=fit$cell_labels),con,row.names=FALSE);close(con)
write_json(list(status='completed',method='SCINA',n_cells=n,library='CM2_glioma_other',rm_overlap=1L,
  solver=solver,native_error=native_error,n_supported_panels=length(sig),
  normalization='Seurat LogNormalize 10000; no geometry required by SCINA',
  predictions_sha256=digest(file=file.path(dest,'predictions.csv.gz'),algo='sha256'),
  source_sha256=digest(file=script,algo='sha256'),job=Sys.getenv('SLURM_JOB_ID')),
  file.path(dest,'fit_manifest.json'),pretty=TRUE,auto_unbox=TRUE)
