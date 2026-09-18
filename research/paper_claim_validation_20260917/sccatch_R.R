# Official scCATCH3.2.2. Union-gene DEG reuse is independently checked against
# library-specific native findmarkergene on the smallest real pilot.
stopifnot(nzchar(Sys.getenv('SLURM_JOB_ID')))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(digest))
root<-'/fs/scratch/PCON0080/yimin/dgscrna'
base<-file.path(root,'results/hvg_ptc_20260916_v1/paper_claim_validation_20260917')
.libPaths(c(file.path(base,'vendor_R'),.libPaths()))
suppressPackageStartupMessages(library(scCATCH))
a<-commandArgs(trailingOnly=TRUE);sample<-a[[1]];mode<-a[[2]]
budget<-if(length(a)>2L)a[[3]] else 'hvg2000'
route<-if(length(a)>3L)a[[4]] else 'PCA30_SNN'
dest<-file.path(base,'comparators/scCATCH',sample,budget,route)
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
if(file.exists(file.path(dest,paste0(toupper(mode),'_COMPLETE'))))quit(status=0)
obj<-readRDS(file.path(base,'GBM',sample,'hvg2000/expression_PCA30.rds'))
X<-obj[['RNA']]@data
cl<-read.csv(file.path(base,'GBM',sample,budget,route,'clusters.csv'),colClasses='character')
stopifnot(identical(cl$cell_id,colnames(X)))
roster<-fromJSON(file.path(base,'markers/scCATCH/manifest.json'),simplifyVector=FALSE)
markers<-lapply(roster$libraries,function(e)read.csv(gzfile(e$path),stringsAsFactors=FALSE))
names(markers)<-vapply(roster$libraries,function(e)e$library,character(1))
union_genes<-unique(unlist(lapply(markers,function(m)m$gene),use.names=FALSE))
X<-X[intersect(rownames(X),union_genes),,drop=FALSE]
cat('SCCATCH_INPUT',sample,budget,route,nrow(X),ncol(X),'\n');flush.console()
native_object<-function(marker,pvalue=.1) {
  ob<-createscCATCH(X,as.character(cl$cluster))
  findmarkergene(ob,if_use_custom_marker=TRUE,marker=marker,use_method='1',
    cell_min_pct=.25,logfc=.25,pvalue=pvalue,verbose=FALSE)
}
pred<-data.frame(cell_id=colnames(X),stringsAsFactors=FALSE);arms<-list();checks<-list()
started<-Sys.time();nclusters<-length(unique(cl$cluster))
union_marker<-data.frame(gene=union_genes,celltype='shared_universe',pmid='union_for_DEG_only',
  subtype1=NA_character_,subtype2=NA_character_,subtype3=NA_character_)
all_deg<-data.frame()
if(nclusters>1L) {
  union_obj<-native_object(union_marker,.1)
  all_deg<-union_obj@markergene
  saveRDS(all_deg,file.path(dest,'official_union_DEG.rds'),compress=FALSE)
}
canonical<-function(d) {
  if(!nrow(d))return(character())
  cols<-c('cluster','gene','comp_cluster','pct','logfc','pvalue')
  d<-d[,cols,drop=FALSE];rownames(d)<-NULL
  sort(apply(d,1,paste,collapse='|'))
}
thresholds<-c(default=.05,permissive=.1,strict=.01)
for(i in seq_along(markers)) {
  marker<-markers[[i]];library<-names(markers)[i]
  marker<-marker[marker$gene %in% rownames(X),,drop=FALSE]
  for(cut in names(thresholds)) {
    pvalue<-thresholds[[cut]];d<-all_deg
    if(nrow(d)) {
      d<-d[d$gene %in% marker$gene & d$pvalue<pvalue,,drop=FALSE]
      if(nrow(d)) {
        key<-paste(d$cluster,d$gene,sep='\034');counts<-table(key)
        d<-d[key %in% names(counts[counts>=nclusters-1L]),,drop=FALSE]
      }
    }
    aid<-paste(sprintf('L%02d',i-1L),cut,sep='__');status<-'completed'
    calls<-rep('Unknown',ncol(X))
    if(nclusters<2L)status<-'no_contrasting_clusters'
    else if(length(unique(marker$gene))<2L)status<-'no_supported_marker_universe'
    else if(!nrow(d))status<-'no_significant_marker_genes'
    if(mode=='audit' && nclusters>1L && length(unique(marker$gene))>=2L) {
      direct<-tryCatch(native_object(marker,pvalue),error=function(e)e)
      if(inherits(direct,'error')) {
        # The official function's nrow(NULL) branch fails if no gene passes.
        # Treat only that exact empty-result case as a no-marker status.
        stopifnot(nrow(d)==0,grepl('argument is of length zero',conditionMessage(direct),fixed=TRUE))
        checks[[aid]]<-list(exact_DEG=TRUE,official_empty_DEG_error=conditionMessage(direct))
      } else {
        stopifnot(identical(canonical(d),canonical(direct@markergene)))
        checks[[aid]]<-list(exact_DEG=TRUE)
      }
    }
    if(nrow(d)>0L) {
      ob<-createscCATCH(X,as.character(cl$cluster));ob@markergene<-d;ob@marker<-marker
      ob<-findcelltype(ob,verbose=FALSE)
      annotation<-ob@celltype
      assigned<-as.character(annotation$cell_type[match(cl$cluster,annotation$cluster)])
      assigned[is.na(assigned)|assigned=='NA'|assigned=='']<-'Unknown'
      calls<-assigned
      if(mode=='audit' && !inherits(direct,'error')) {
        native<-findcelltype(direct,verbose=FALSE)
        other<-as.character(native@celltype$cell_type[match(cl$cluster,native@celltype$cluster)])
        other[is.na(other)|other=='NA'|other=='']<-'Unknown'
        stopifnot(identical(calls,other));checks[[aid]]$exact_final_calls<-TRUE
      }
    }
    pred[[aid]]<-calls
    arms[[aid]]<-list(method='scCATCH',library=library,budget=budget,route=route,cutoff=cut,
      status=status,pvalue=pvalue,use_method='1',cell_min_pct=.25,logfc=.25,
      cluster_sha256=digest(file=file.path(base,'GBM',sample,budget,route,'clusters.csv'),algo='sha256'))
  }
  cat('SCCATCH_LIBRARY',sample,library,mode,'\n');flush.console()
}
con<-gzfile(file.path(dest,paste0(mode,'_predictions.csv.gz')),'wt');write.csv(pred,con,row.names=FALSE);close(con)
m<-list(status='completed',method='scCATCH',sample=sample,mode=mode,arms=arms,checks=checks,
  version=as.character(packageVersion('scCATCH')),n_cells=ncol(X),thresholds=thresholds,
  input='full log1p library-size normalized RNA; matched marker genes; native clusters including the HDBSCAN noise partition',
  ambiguous_calls='All tied native labels retained as returned by the official function; mapping must preserve ambiguity.',
  elapsed_seconds=as.numeric(difftime(Sys.time(),started,units='secs')),job=Sys.getenv('SLURM_JOB_ID'),
  predictions_sha256=digest(file=file.path(dest,paste0(mode,'_predictions.csv.gz')),algo='sha256'),
  source_sha256=digest(file=sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1]),algo='sha256'))
write_json(m,file.path(dest,paste0(mode,'_manifest.json')),pretty=TRUE,auto_unbox=TRUE)
writeLines(digest(file=file.path(dest,paste0(mode,'_manifest.json')),algo='sha256'),file.path(dest,paste0(toupper(mode),'_COMPLETE')))
