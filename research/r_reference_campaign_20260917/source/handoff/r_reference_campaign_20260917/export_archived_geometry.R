stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
suppressPackageStartupMessages(library(Seurat))
root<-'/fs/scratch/PCON0080/yimin/dgscrna'
base<-file.path(root,'results/hvg_ptc_20260916_v1')
out<-file.path(base,'r_reference_campaign_20260917/PTC_archived_CCA2000')
e<-new.env();load(file.path(base,'ptc_recovery/archive/tcr/rawdata/integrated_data_final_annotation.Rdata'),envir=e)
obj<-e$integrated_data
cells<-read.csv(file.path(out,'seurat_clusters/cells.csv'),stringsAsFactors=FALSE)
stopifnot(identical(cells$original_cell_id,colnames(obj)))
for(red in c('pca','umap')) {
 z<-Embeddings(obj,red);stopifnot(identical(rownames(z),colnames(obj)))
 rownames(z)<-cells$cell_id
 write.csv(z,file.path(out,if(red=='pca')'PCA30.csv' else 'UMAP2.csv'))
}
writeLines('original saved embeddings; no recomputation',file.path(out,'GEOMETRY_EXPORTED'))
