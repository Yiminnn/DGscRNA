#!/usr/bin/env Rscript
# Independently replay historical terminal annotations in R, then export fixed-QC counts.
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(jsonlite))
root <- '/fs/scratch/PCON0080/yimin/dgscrna'
base <- file.path(root,'results/hvg_ptc_20260916_v1/ptc_recovery')
dest <- file.path(base,'R_baseline_replay')
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
stopifnot(file.exists(file.path(base,'table_reconciliation/COMPLETE')))
report <- list(status='running',job=Sys.getenv('SLURM_JOB_ID'),R=R.version.string)
write_report <- function() jsonlite::write_json(report,file.path(dest,'replay.json'),pretty=TRUE,auto_unbox=TRUE,na='null')
write_report()
source_path <- file.path(base,'archive/tcr/rawdata/integrated_data_final_annotation.Rdata')
historical <- new.env(parent=emptyenv())
loaded <- load(source_path,envir=historical)
stopifnot(identical(loaded,'integrated_data'))
obj <- historical$integrated_data
meta <- methods::slot(obj,'meta.data')
original <- read.delim(file.path(base,'archive/tcr/rawdata/metadata.txt'),check.names=FALSE)
sample_map <- setNames(original$Sample,original$sc_ID)
patient_map <- setNames(original$Patient,original$Sample)
canonical <- paste0(sample_map[meta$sample.name],'_',sub('-.*','',rownames(meta)))
stopifnot(!anyDuplicated(canonical),length(canonical)==92404L,!anyNA(canonical))
s2 <- read.csv(gzfile(file.path(base,'table_reconciliation/S2_as_supplied.csv.gz')),row.names=1,check.names=FALSE)
stopifnot(setequal(canonical,rownames(s2)))
s2 <- s2[canonical,,drop=FALSE]
correspondence <- c('Cell_type_annotation_CellMarker2.0'='cell_type_annotation_CellMarker2.0',
                   'DG_scRNA_Finalized_Cell_Types'='DGCyTOF_Finalized_Cell_Types',
                   'Cell_types_general'='cell_types_general')
agreement <- lapply(names(correspondence),function(target) {
  source <- unname(correspondence[[target]])
  equal <- as.character(s2[[target]])==as.character(meta[[source]])
  stopifnot(!anyNA(equal))
  data.frame(target=target,source=source,n=length(equal),n_equal=sum(equal),n_mismatch=sum(!equal))
})
agreement <- do.call(rbind,agreement)
write.csv(agreement,file.path(dest,'S2_R_independent_column_agreement.csv'),row.names=FALSE)
stopifnot(all(agreement$n_mismatch==0))
generalize_legacy <- function(x) {
  if (x %in% c('Undecided','Unknown','No_Annotation')) return(x)
  fields <- strsplit(x,'+',fixed=TRUE)[[1]]
  discard <- if (fields[1]=='NCOMMREFF') 1L else if (fields[1]=='cancer') 3L else 2L
  stopifnot(length(fields)>discard)
  paste(fields[seq.int(discard+1L,length(fields))],collapse='+')
}
native <- as.character(meta$DGCyTOF_Finalized_Cell_Types)
dictionary <- setNames(vapply(unique(native),generalize_legacy,character(1)),unique(native))
general <- unname(dictionary[native])
stopifnot(identical(general,as.character(meta$cell_types_general)))
write.csv(data.frame(native=names(dictionary),general=unname(dictionary)),
          file.path(dest,'native_to_general_exact_mapping.csv'),row.names=FALSE)
initial <- as.character(meta$cell_type_annotation_CellMarker2.0)
initial_unknown <- initial=='Undecided'
report$annotation_replay <- list(n_cells=length(native),n_samples=length(unique(meta$sample.name)),
  n_patients=length(unique(meta$patient.name)),n_S2_native_mismatches=0L,n_S2_general_mismatches=0L,
  n_initial_undecided=sum(initial_unknown),n_changed_previously_known=sum(!initial_unknown & initial!=native),
  n_terminal_unknown=sum(native %in% c('Undecided','Unknown','No_Annotation')),
  source='archived terminal predictions; no historical weight checkpoint recovered',
  stochastic_retraining_exactness='not claimed')
report$historical_geometry <- list(seurat_version=as.character(methods::slot(obj,'version')),
  integrated_n_features=nrow(obj[['integrated']]),n_cells=ncol(obj),
  n_samples_in_CCA=length(methods::slot(methods::slot(obj,'commands')$FindIntegrationAnchors,'assay.used')),
  PCA_n_components=ncol(Embeddings(obj,'pca')),UMAP_n_components=ncol(Embeddings(obj,'umap')),
  CCA_grouping='all eight samples in the saved S2 checkpoint; new requested grouping is an explicit intervention')
report$status <- 'archived_terminal_annotation_replay_verified'
write_report()
writeLines(capture.output(sessionInfo()),file.path(dest,'sessionInfo.txt'))

# Counts are exported only after the R-side per-cell replay has passed.
counts <- methods::slot(obj[['RNA']],'counts')
stopifnot(identical(colnames(counts),rownames(meta)))
stopifnot(all(is.finite(counts@x)),all(counts@x>=0),all(counts@x==round(counts@x)))
report$raw_counts <- list(n_features=nrow(counts),n_cells=ncol(counts),
  n_stored_values=length(counts@x),finite_nonnegative_integer=TRUE,
  cell_selection='exact historical 92404 QC-retained cells; no new doublet exclusion',
  gene_selection='all RNA assay features; feature filtering is deferred to explicit analysis stages')
groups <- list(MTN=c('MT-1','MT-2','N-1','N-2'),TUT=c('TU-1','TU-2','T-1','T-2'))
input_dir <- file.path(base,'fixed_QC_inputs');dir.create(input_dir,showWarnings=FALSE)
qc <- meta[,c('nCount_RNA','nFeature_RNA','percent.mt','sample.name','patient.name','group.name'),drop=FALSE]
qc$sample_id <- unname(sample_map[qc$sample.name])
qc$patient_id <- unname(patient_map[qc$sample_id])
rownames(qc) <- canonical
colnames(counts) <- canonical
write.csv(qc,file.path(input_dir,'cells_QC_sample_patient.csv'))
for (name in names(groups)) {
  cat('Exporting fixed-QC count input',name,'\n');flush.console()
  keep <- which(qc$sample_id %in% groups[[name]])
  payload <- list(counts=counts[,keep,drop=FALSE],qc_metadata=qc[keep,,drop=FALSE],
                  historical_source=source_path,analysis_group=name,
                  annotation_columns_removed=TRUE)
  saveRDS(payload,file.path(input_dir,paste0(name,'.counts.rds')),compress=FALSE)
  write.csv(data.frame(cell_id=colnames(payload$counts)),file.path(input_dir,paste0(name,'.cells.csv')),row.names=FALSE)
  report$exported_groups[[name]] <- list(n_cells=length(keep),samples=groups[[name]],n_genes=nrow(counts))
  rm(payload);gc()
}
write_report()
writeLines(format(Sys.time(),tz='UTC'),file.path(dest,'COMPLETE'))
cat('R per-cell terminal annotation replay passed and fixed-QC inputs exported.\n')
