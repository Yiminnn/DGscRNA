#!/usr/bin/env Rscript
# Reuse the frozen, independently verified scorer; dispatch only a new analysis scope.
source_path <- Sys.getenv('PTC_SCORE_SOURCE')
parsed <- parse(source_path)
stopifnot(is.call(parsed[[length(parsed)-1L]]),identical(parsed[[length(parsed)-1L]][[1]],as.name('if')))
stopifnot(grepl('mode %in%',paste(deparse(parsed[[length(parsed)-1L]][[2]]),collapse=' '),fixed=TRUE))
for(expression in parsed[seq_len(length(parsed)-2L)]) eval(expression)
sample <- c('MT-1','MT-2','N-1','N-2','TU-1','TU-2','T-1','T-2')[[task_id+1L]]
group <- if(task_id<4L)'MTN' else 'TUT'
obj <- readRDS(file.path(base,'samples',sample,'RNA_normalized_VST.rds'))
gd <- file.path(base,'matched_samples',sample)
stopifnot(file.exists(file.path(gd,'GEOMETRY_COMPLETE')))
for(space in c('PCA30','UMAP2')) for(method in c('SNN','HDBSCAN_R')) {
  context <- list(sample=sample,group=group,correction='per_sample_NONE',space=space,
      clusterer=method,family='matched_R_single',seed=42,
      dispatch_sha256=digest(file=Sys.getenv('PTC_MATCHED_SCORE'),algo='sha256'))
  score_one(obj,file.path(gd,paste0(space,'_',method)),'RNA',names(libraries),c('none','mean','0.5'),context)
}
cat('Matched single-sample scoring complete\n')
