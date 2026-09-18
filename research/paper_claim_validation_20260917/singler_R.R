# SingleR2.8.0 with labelled training patients; exact nearest-neighbor search,
# built-in reference aggregation and fine tuning. Test labels are not loaded.
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
suppressPackageStartupMessages(library(Matrix));suppressPackageStartupMessages(library(SingleR))
suppressPackageStartupMessages(library(jsonlite));suppressPackageStartupMessages(library(digest))
suppressPackageStartupMessages(library(BiocParallel))
root<-'/fs/scratch/PCON0080/yimin/dgscrna';base<-file.path(root,'results/hvg_ptc_20260916_v1/paper_claim_validation_20260917')
args<-commandArgs(trailingOnly=TRUE);fold<-as.integer(args[[1]])
rd<-file.path(base,'reference_inputs',paste0('fold',fold));m<-fromJSON(file.path(rd,'manifest.json'))
stopifnot(readLines(file.path(rd,'COMPLETE'))==digest(file=file.path(rd,'manifest.json'),algo='sha256'))
load_binary<-function(path,ncells,ngenes,nnz,cells,genes) {
  read<-function(n,what,size,count) {con<-file(file.path(path,n),'rb');on.exit(close(con));readBin(con,what,n=count,size=size,endian='little')}
  new('dgCMatrix',x=read('x.bin','numeric',8,nnz),i=read('i.bin','integer',4,nnz),
    p=read('p.bin','integer',4,ncells+1L),Dim=as.integer(c(ngenes,ncells)),Dimnames=list(genes,cells))
}
robs<-read.csv(file.path(rd,'reference_cells.csv'),stringsAsFactors=FALSE);genes<-read.csv(file.path(rd,'genes.csv'))$gene
ref<-load_binary(rd,m$n_cells,m$n_genes,m$nnz,robs$reference_id,genes)
tests<-read.csv(file.path(rd,'test_samples.csv'),stringsAsFactors=FALSE)
stopifnot(!any(tests$patient %in% robs$patient))
BPPARAM<-MulticoreParam(workers=min(4L,as.integer(Sys.getenv('SLURM_CPUS_PER_TASK','4'))),RNGseed=42L)
set.seed(42)
# Use only genes available in every target sample in this fold. Gene identities
# are label-free input metadata; expression/labels never enter reference fitting.
common<-genes
for(sample in tests$sample)common<-intersect(common,read.csv(file.path(base,'inputs',sample,'genes.csv'))$gene)
stopifnot(length(common)>1000L)
ref<-ref[common,,drop=FALSE]
trained<-trainSingleR(ref,labels=robs$L1,test.genes=common,aggr.ref=TRUE,
  de.method='classic',BPPARAM=BPPARAM)
cat('SINGLER_TRAINED',fold,ncol(ref),nrow(ref),'\n');flush.console()
outfold<-file.path(base,'comparators/SingleR',paste0('fold',fold));dir.create(outfold,recursive=TRUE,showWarnings=FALSE)
writeLines(common,file.path(outfold,'reference_gene_universe.txt'))
for(j in seq_len(nrow(tests))) {
  sample<-tests$sample[j];dest<-file.path(base,'comparators/SingleR',sample)
  dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  if(file.exists(file.path(dest,'COMPLETE')))next
  src<-file.path(base,'inputs',sample);im<-fromJSON(file.path(src,'input_manifest.json'))
  stopifnot(!(tests$patient[j] %in% m$training_patients))
  cells<-read.csv(file.path(src,'cells_fit.csv'))$cell_id;gg<-read.csv(file.path(src,'genes.csv'))$gene
  X<-load_binary(src,im$n_cells,im$n_genes,im$nnz,cells,gg)
  totals<-Matrix::colSums(X);stopifnot(all(totals>0))
  X@x<-log1p(X@x*rep(10000/totals,diff(X@p)));X<-X[common,,drop=FALSE]
  start<-Sys.time();fit<-classifySingleR(X,trained,fine.tune=TRUE,prune=TRUE,BPPARAM=BPPARAM)
  raw<-as.character(fit$labels);pruned<-as.character(fit$pruned.labels);pruned[is.na(pruned)]<-'Unknown'
  stopifnot(identical(rownames(fit),cells),all(raw %in% unique(robs$L1)))
  pred<-data.frame(cell_id=cells,default=pruned,unpruned=raw,stringsAsFactors=FALSE)
  con<-gzfile(file.path(dest,'predictions.csv.gz'),'wt');write.csv(pred,con,row.names=FALSE);close(con)
  saveRDS(fit,file.path(dest,'scores.rds'),compress=FALSE)
  report<-list(status='completed',method='SingleR',version=as.character(packageVersion('SingleR')),sample=sample,
    fold=fold,n_cells=ncol(X),n_reference_cells=ncol(ref),n_common_genes=nrow(X),reference_labels=sort(unique(robs$L1)),
    reference_manifest_sha256=digest(file=file.path(rd,'manifest.json'),algo='sha256'),
    no_heldout_patient_in_reference=TRUE,supervision='Curated training-patient labels, primary97 reference only',
    primary='default pruned per-cell calls; all test cells retained; unpruned calls reported as sensitivity',
    aggr_ref=TRUE,de_method='classic',fine_tune=TRUE,seed=42L,
    elapsed_seconds=as.numeric(difftime(Sys.time(),start,units='secs')),job=Sys.getenv('SLURM_JOB_ID'),
    predictions_sha256=digest(file=file.path(dest,'predictions.csv.gz'),algo='sha256'),
    source_sha256=digest(file=sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]),algo='sha256'))
  write_json(report,file.path(dest,'manifest.json'),pretty=TRUE,auto_unbox=TRUE)
  writeLines(digest(file=file.path(dest,'manifest.json'),algo='sha256'),file.path(dest,'COMPLETE'))
  cat('SINGLER_PREDICTED',sample,ncol(X),'\n');flush.console()
}
write_json(list(status='completed',fold=fold,samples=tests$sample,job=Sys.getenv('SLURM_JOB_ID')),
  file.path(outfold,'manifest.json'),pretty=TRUE,auto_unbox=TRUE)
