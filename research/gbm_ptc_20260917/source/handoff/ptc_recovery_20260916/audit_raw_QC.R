#!/usr/bin/env Rscript
# QC/count provenance audit only. Never change the fixed matched-experiment cell set.
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest))
root <- '/fs/scratch/PCON0080/yimin/dgscrna'
base <- file.path(root,'results/hvg_ptc_20260916_v1/ptc_experiments')
recovery <- file.path(root,'results/hvg_ptc_20260916_v1/ptc_recovery')
samples <- c('MT-1','MT-2','N-1','N-2','TU-1','TU-2','T-1','T-2')
i <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'));sample <- samples[[i+1L]]
group <- if(i<4)'MTN' else 'TUT'
meta <- read.delim(file.path(recovery,'archive/tcr/rawdata/metadata.txt'))
scid <- meta$sc_ID[meta$Sample==sample];stopifnot(length(scid)==1L)
dir <- file.path(root,'data_ptc/tcr/rawdata/GEX',scid)
raw <- Read10X(dir,gene.column=2,unique.features=TRUE)
if(is.list(raw))raw <- raw[['Gene Expression']]
stopifnot(inherits(raw,'sparseMatrix'),all(raw@x>=0),all(raw@x==round(raw@x)))
raw_gene_names <- rownames(raw)
rownames(raw) <- gsub('_','-',raw_gene_names,fixed=TRUE)
stopifnot(!anyDuplicated(rownames(raw)))
gene_name_adapter <- data.frame(raw_symbol=raw_gene_names,seurat_symbol=rownames(raw))
gene_name_adapter <- gene_name_adapter[gene_name_adapter$raw_symbol!=gene_name_adapter$seurat_symbol,,drop=FALSE]
colnames(raw) <- paste0(sample,'_',sub('-.*','',colnames(raw)))
stopifnot(!anyDuplicated(colnames(raw)))
payload <- readRDS(file.path(recovery,'fixed_QC_inputs',paste0(group,'.counts.rds')))
selected <- payload$qc_metadata$sample_id==sample
saved <- payload$counts[,selected,drop=FALSE];qc <- payload$qc_metadata[selected,,drop=FALSE]
stopifnot(all(colnames(saved) %in% colnames(raw)))
rawgenes <- rownames(raw)
gene3 <- Matrix::rowSums(raw>0)>=3L
original_gene3 <- raw[gene3,,drop=FALSE]
summarize <- function(counts) {
 n <- Matrix::colSums(counts);f <- Matrix::colSums(counts>0)
 mt <- Matrix::colSums(counts[grepl('^MT-',rownames(counts)),,drop=FALSE])
 data.frame(cell_id=colnames(counts),nCount=n,nFeature=f,percent_mt=100*mt/pmax(n,1),row.names=NULL)
}
raw_qc <- summarize(raw);filtered_qc <- summarize(original_gene3)
stopifnot(identical(raw_qc$cell_id,filtered_qc$cell_id))
records <- data.frame(cell_id=raw_qc$cell_id,sample=sample,
  retained_in_historical_S2=raw_qc$cell_id %in% colnames(saved),
  raw_nCount=raw_qc$nCount,raw_nFeature=raw_qc$nFeature,raw_percent_mt=raw_qc$percent_mt,
  gene3_nCount=filtered_qc$nCount,gene3_nFeature=filtered_qc$nFeature,gene3_percent_mt=filtered_qc$percent_mt)
records$reference_constructor_pass <- records$gene3_nFeature>=200
records$reference_QC_pass <- records$reference_constructor_pass & records$gene3_nFeature>200 & records$gene3_percent_mt<15
records$raw_reference_QC_pass <- records$raw_nFeature>200 & records$raw_percent_mt<15
common <- intersect(rownames(saved),rownames(raw))
stopifnot(length(common)>30000L)
delta <- raw[common,colnames(saved),drop=FALSE]-saved[common,,drop=FALSE]
comparison <- data.frame(cell_id=colnames(saved),
   changed_genes_raw_vs_archived=Matrix::colSums(delta!=0),
   absolute_count_difference=Matrix::colSums(abs(delta)),
   archived_nCount=qc$nCount_RNA,archived_nFeature=qc$nFeature_RNA,archived_percent_mt=qc$percent.mt)
q <- filtered_qc[match(colnames(saved),filtered_qc$cell_id),,drop=FALSE]
comparison$gene3_nCount_minus_archived <- q$nCount-qc$nCount_RNA
comparison$gene3_nFeature_minus_archived <- q$nFeature-qc$nFeature_RNA
comparison$gene3_percent_mt_minus_archived <- q$percent_mt-qc$percent.mt
raw_q <- raw_qc[match(colnames(saved),raw_qc$cell_id),,drop=FALSE]
comparison$raw_nCount_minus_archived <- raw_q$nCount-qc$nCount_RNA
comparison$raw_nFeature_minus_archived <- raw_q$nFeature-qc$nFeature_RNA
comparison$raw_percent_mt_minus_archived <- raw_q$percent_mt-qc$percent.mt
missing_in_raw <- setdiff(rownames(saved),rownames(raw))
missing_in_saved <- setdiff(rownames(raw),rownames(saved))
dest <- file.path(base,'raw_QC_audit_symbol_normalized');dir.create(dest,showWarnings=FALSE)
write.csv(gene_name_adapter,file.path(dest,paste0(sample,'.Seurat_gene_name_adapter.csv')),row.names=FALSE)
con <- gzfile(file.path(dest,paste0(sample,'.all_raw_cells.csv.gz')),'wt');write.csv(records,con,row.names=FALSE);close(con)
con <- gzfile(file.path(dest,paste0(sample,'.retained_count_comparison.csv.gz')),'wt');write.csv(comparison,con,row.names=FALSE);close(con)
nqc <- sum(records$reference_QC_pass)
m <- list(sample=sample,group=group,n_raw_cells=ncol(raw),n_raw_genes=nrow(raw),
 n_genes_detected_in3_cells=sum(gene3),n_after_reference_constructor=sum(records$reference_constructor_pass),
 n_after_reference_QC=nqc,n_archived_S2=ncol(saved),
 n_archived_cells_failing_candidate_QC=sum(records$retained_in_historical_S2 & !records$reference_QC_pass),
 n_candidate_QC_cells_absent_from_S2=sum(records$reference_QC_pass & !records$retained_in_historical_S2),
 candidate_expected_doublets_round075=round(0.075*nqc),
 candidate_expected_after_doublets=nqc-round(0.075*nqc),
 n_after_raw_reference_QC=sum(records$raw_reference_QC_pass),
 raw_candidate_expected_after_doublets=sum(records$raw_reference_QC_pass)-round(0.075*sum(records$raw_reference_QC_pass)),
 n_archived_cells_failing_raw_candidate_QC=sum(records$retained_in_historical_S2 & !records$raw_reference_QC_pass),
 n_retained_cells_raw_count_exact=sum(comparison$changed_genes_raw_vs_archived==0),
 n_retained_cells_raw_nCount_exact=sum(comparison$raw_nCount_minus_archived==0),
 n_retained_cells_raw_nFeature_exact=sum(comparison$raw_nFeature_minus_archived==0),
 n_retained_cells_raw_mt_exact=sum(abs(comparison$raw_percent_mt_minus_archived)<1e-10),
 n_retained_cells_gene3_nCount_exact=sum(comparison$gene3_nCount_minus_archived==0),
 n_retained_cells_gene3_nFeature_exact=sum(comparison$gene3_nFeature_minus_archived==0),
 n_retained_cells_gene3_mt_exact=sum(abs(comparison$gene3_percent_mt_minus_archived)<1e-10),
 archived_genes_absent_from_raw=missing_in_raw,raw_genes_absent_from_archived=missing_in_saved,
 all_archived_genes_compared=length(missing_in_raw)==0L && length(missing_in_saved)==0L,
 n_Seurat_underscore_gene_name_changes=nrow(gene_name_adapter),
 candidate_definition='Archived reference vignette CreateSeuratObject min.cells3/min.features200; source.R nFeature>200 & percent.mt<15; expected DoubletFinder7.5%',
 limitation='Reference scripts include other cohorts. Candidate QC/count compatibility is not proof of original stochastic doublet labels; no original pK/weights recovered. No fixed experimental cell set is changed.',
 files=lapply(c('matrix.mtx.gz','features.tsv.gz','barcodes.tsv.gz'),function(name)
   list(path=file.path(dir,name),sha256=digest(file=file.path(dir,name),algo='sha256'))),
 job=Sys.getenv('SLURM_JOB_ID'))
write_json(m,file.path(dest,paste0(sample,'.json')),pretty=TRUE,auto_unbox=TRUE)
writeLines(digest(file=file.path(dest,paste0(sample,'.json')),algo='sha256'),file.path(dest,paste0(sample,'.COMPLETE')))
cat(sample,'raw count/QC provenance audit complete\n')
