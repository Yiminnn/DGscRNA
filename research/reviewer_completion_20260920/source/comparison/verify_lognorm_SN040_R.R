# Verify actual CSV supplied to scDeepSort against original R RNA assay, every value.
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
suppressPackageStartupMessages({library(Seurat);library(data.table);library(jsonlite);library(digest)})
root<-'/fs/scratch/PCON0080/yimin/dgscrna'
old<-file.path(root,'results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/GBM')
out<-file.path(root,'results/hvg_ptc_20260916_v1/reviewer_completion_20260920/comparison/scDeepSort_LogNormalize')
proof<-list()
for(sample in c('SN040')) {
  path<-file.path(out,sample,'lognorm_model_input.csv')
  df<-fread(path,check.names=FALSE)
  genes<-as.character(df[[1]])
  matrix<-as.matrix(df[,-1,with=FALSE]);rownames(matrix)<-genes
  rm(df);gc()
  obj<-readRDS(file.path(old,sample,'hvg2000/expression_PCA30.rds'))
  stopifnot(identical(colnames(matrix),colnames(obj)),!anyDuplicated(genes),all(genes %in% rownames(obj[['RNA']]@data)))
  reference<-as.matrix(obj[['RNA']]@data[genes,colnames(matrix),drop=FALSE])
  delta<-max(abs(matrix-reference))
  stopifnot(is.finite(delta),delta<1e-12)
  proof[[sample]]<-list(sample=sample,n_cells=ncol(matrix),n_genes=nrow(matrix),
    n_values_checked=length(matrix),max_absolute_difference=delta,passed=TRUE,
    actual_predictor_csv_sha256=digest(file=path,algo='sha256'),
    R_expression_object_sha256=digest(file=file.path(old,sample,'hvg2000/expression_PCA30.rds'),algo='sha256'))
  cat('OFFICIAL_INPUT_FULL_R_PARITY',sample,ncol(matrix),nrow(matrix),delta,'\n')
  rm(obj,reference,matrix);gc()
}
write_json(list(status='passed',job=Sys.getenv('SLURM_JOB_ID'),samples=proof,
  definition='Actual real-valued predictor CSV versus native R RNA data assay, every shared gene and cell',
  source_sha256=digest(file=sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]),algo='sha256')),
  file.path(out,'R_full_input_parity_SN040.json'),pretty=TRUE,auto_unbox=TRUE)
