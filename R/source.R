BioCManagerLibraries <- c()
# requiredLibraries <- c(BioCManagerLibraries, 
#                        "Seurat", 
#                        "SeuratData", 
#                        "future", 
#                        #"ggpubr", 
#                        'gridExtra', 
#                        "tidyverse", 
#                        "DoubletFinder", 
#                        "SeuratDisk", 
#                        "Rcpp",
#                        "dbscan",
#                        "ggplot2",
#                        "openxlsx",
#                        'stringr',
#                        'nloptr',
#                        'monocle3',
#                        'scCCESS')
# 
# #remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
# #remotes::install_github("mojaveazure/seurat-disk")
# #devtools::install_github('cole-trapnell-lab/monocle3')
# #devtools::install_github('PYangLab/scCCESS')
# #install.packages('clusterCrit_1.2.8.tar.gz', repos = NULL, type ='source')
# for (packageName in requiredLibraries){
#   if (!is.element(packageName, installed.packages()[,1])){
#     print(paste("Installing package: ", packageName))
#     if (packageName %in% BioCManagerLibraries) {
#       BiocManager::install(packageName, INSTALL_opts = '--no-lock') # type="binary"
#     }
#     else {
#       install.packages(packageName, dependencies = TRUE, INSTALL_opts = '--no-lock') #type="binary")#
#       #install.packages("table1", dependencies = TRUE, INSTALL_opts = '--no-lock')
# 
#     }
# 
#   }
# 
#   suppressMessages(library(packageName, character.only = TRUE))
#   print(paste("Loaded package: ", packageName))
# 
# }

requiredLibraries <- c("Seurat", 
                       "SeuratData", 
                       "future", 
                       'gridExtra', 
                       "tidyverse", 
                       "DoubletFinder", 
                       "SeuratDisk", 
                       "Rcpp",
                       "dbscan",
                       "ggplot2",
                       "openxlsx",
                       'stringr',
                       'nloptr',
                       'monocle3',
                       'scCCESS',
                       'reticulate',
                       'SeuratWrappers',
                       'R.utils',
                       'tensorflow',
                       'clue',
                       "dplyr",
                       "HGNChelper",
                       'pROC',
                       'MLmetrics')

for (packageName in requiredLibraries){
  suppressMessages(library(packageName, character.only = TRUE))
}


# library('SeuratData')
# library('Seurat')
# InstallData("pbmc3k")

# data("pbmc3k")
# seurat_object@project.name

#' Preprocessing for a single cell Seurat object
#'
#' @param seurat_object A single cell Seurat Object
#' @return Returns a list:
#'
#' * `data` - preprocessed single cell data
#' * `qc_info` - Info regarding how many cells were filtered and at what stage
#' * `qc_plot` - Initial qc plot of gene expression in cells
preprocessing <- function(seurat_object){
  # preprocess the input seurat object
  print(paste("Total Number of samples:", ncol(seurat_object)))
  all_available_cells <- ncol(seurat_object)
  
  seurat_object[["percent.mt"]]<-PercentageFeatureSet(seurat_object,pattern="^MT-")
  
  # Get the plot for QC
  qc_plot <- VlnPlot(seurat_object,features=c("nCount_RNA", "nFeature_RNA","percent.mt"),ncol=3)
  
  # Apple QC, get subset of the features
  seurat_object <-subset(seurat_object, subset = nFeature_RNA > 200 & percent.mt < 15)
  
  after_qc_cells <- ncol(seurat_object)
  
  print(paste("Total Number of samples (after QC):", ncol(seurat_object)))
  
  ## Second, Normalize data and calculate variable features using VST
  seurat_object <- NormalizeData(seurat_object, 
                                 normalization.method = "LogNormalize",
                                 scale.factor = 10000)
  
  # Calculate variable features
  seurat_object <- FindVariableFeatures(seurat_object, 
                                        selection.method = "vst",
                                        nfeatures = 2000)
  # Scale Data and run PCA
  seurat_object <- ScaleData(seurat_object, vars.to.regress = c("percent.mt"))
  
  seurat_object <- RunPCA(seurat_object, npcs = 30, verbose = FALSE)
  
  # add meta data, group name, sample name, and patient name
  # res <- AddMetaData(object = res,
  #                    metadata = rep(all_data_samples[[i]], length(colnames(res))),
  #                    col.name = "sample.name")
  
  # res <- AddMetaData(object = res,
  #                    metadata = rep(group_data[[i]], length(colnames(res))),
  #                    col.name = "group.name")
  
  # res <- AddMetaData(object = res,
  #                    metadata = rep(patient_data[[i]], length(colnames(res))),
  #                    col.name = "patient.name")
  
  
  
  # doplet removal
  
  seurat_object.sweep <- paramSweep_v3(seurat_object, PCs = 1:10, sct = FALSE)
  seurat_object.stats <- summarizeSweep(seurat_object.sweep, GT = FALSE)
  seurat_object.PK <- find.pK(seurat_object.stats)
  
  best_pK = which(seurat_object.PK$BCmetric == max(seurat_object.PK$BCmetric))
  best_pK <- as.numeric(levels(seurat_object.PK$pK)[best_pK])
  
  nExp_poi <- round(0.075*nrow(seurat_object@meta.data)) # assume 7.5% droplet chance, according to the pipeline
  
  # find droplets and label them as such
  seurat_object <- doubletFinder_v3(seurat_object, PCs = 1:10, pN = 0.25, pK = best_pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  
  # according to classification, remove the droplets as QC once again
  fetched <- FetchData(object = seurat_object, vars = tail(colnames(seurat_object@meta.data), 1))
  seurat_object <- seurat_object[, which(fetched == "Singlet")]
  
  print(paste("Total Number of samples (after QC & doublet removal):", ncol(seurat_object)))
  
  removed_qc_cells <- data.frame("sample" = seurat_object@project.name,
                                 "total_num_cells" = all_available_cells, 
                                 "total_num_cells_after_qc" = after_qc_cells,
                                 "chosen_pk_val" = best_pK,
                                 "total_num_cells_after_qc_droplet" = ncol(seurat_object), 
                                 "total_num_cells_filtered" = all_available_cells - ncol(seurat_object))
  # removed_qc_cells <- rbind(removed_qc_cells, c(seurat_object@project.name,
  #                                               all_available_cells,
  #                                               after_qc_cells,
  #                                               best_pK,
  #                                               ncol(seurat_object),
  #                                               all_available_cells - ncol(seurat_object)))
  # 
  # colnames(removed_qc_cells) <- c("sample", 
  #                                 "total_num_cells", 
  #                                 "total_num_cells_after_qc",
  #                                 "chosen_pk_val",
  #                                 "total_num_cells_after_qc_droplet", 
  #                                 "total_num_cells_filtered")
  
  return(list('data' = seurat_object,
              'qc_info' = removed_qc_cells,
              'qc_plot' = qc_plot))
  
}

#' Preprocessing for a list of single cell seurat objects
#'
#' @param list_of_seurat A list of seurat objects
#' @return Returns a list:
#'
#' * `integrated_data` - preprocessed and batch effect removed data
#' * `list_of_seurat` - Individual seurat information, from preprocessing function
list_preprocessing <- function(list_of_seurat){
  res <- list('data' = c(),
              'qc_info' = c(),
              'qc_plot' = c())
  for (item in panc8){
    # run preprocessing for the 
    object_res <- preprocessing(item)
    res[['data']] <- c(res[['data']], object_res[['data']])
    res[['qc_info']] <- c(res[['qc_info']], object_res[['qc_info']])
    res[['qc_plot']] <- c(res[['qc_plot']], object_res[['qc_plot']])
  }
  
  # now we proceed to run seurat batch effect removal
  print('Starting integration...')
  
  # proceed to integrate them together
  features <- SelectIntegrationFeatures(object.list = res[['data']])
  
  integration_anchors <- FindIntegrationAnchors(object.list = res[['data']], anchor.features = features)
  
  integrated_data <- IntegrateData(anchorset = integration_anchors)
  
  return(list('integrated_data' = integrated_data, 
              'list_of_seurat' = res))
}

#' Apply dimensionality reductions to the integrated seurat object, currently supports PCA and UMAP
#'
#' @param integrated_data A seurat object
#' @return a seurat object, dimensionality reductions
dim_reductions <- function(integrated_data){

  ### proceed with Visualization
  DefaultAssay(integrated_data) <- "integrated"
  
  # Run the standard workflow for visualization and clustering
  integrated_data <- ScaleData(integrated_data, verbose = FALSE)
  integrated_data <- RunPCA(integrated_data, npcs = 30, verbose = FALSE)
  integrated_data <- RunUMAP(integrated_data, reduction = "pca", dims = 1:30)
  
  return(integrated_data)
}

#' Apply clustering strategies to the integrated data with dimensionality reduction. Currently supports Seurat, HDBSCAN, monocle3, and scCCESS
#'
#' @param integrated_data A seurat object
#' @return a seurat object
clustering_strategies <- function(integrated_data, strategy = c('Seurat', 'HDBSCAN', 'monocle3', 'scCCESS')){
  
  for (select_strategy in c(strategy)){
    print(paste('Attempting to run',select_strategy,'...'))
    if(select_strategy == 'Seurat'){
      
      # seurat - PCA
      print('PCA reduction clusters...')
      integrated_data <- FindNeighbors(integrated_data, reduction = "pca", dims = 1:30)
      integrated_data <- FindClusters(integrated_data, resolution = 0.5)

      # seurat - UMAP
      print('UMAP reduction clusters...')
      temp <- FindNeighbors(integrated_data, reduction = "umap", dims = 1:2)
      temp <- FindClusters(temp, resolution = 0.5)
      
      integrated_data <- AddMetaData(integrated_data,
                                     metadata = temp$integrated_snn_res.0.5,
                                     col.name = 'seurat.UMAP_clusters')
      # release object
      rm(temp)
      
    }else if(select_strategy == 'HDBSCAN'){
      
      # hdbscan - PCA
      print('PCA reduction clusters...')
      hdbscan_res <- dbscan::hdbscan(integrated_data@reductions$pca@cell.embeddings,
                                     minPts = 50)
      integrated_data <- AddMetaData(integrated_data,
                                     metadata = hdbscan_res$cluster,
                                     col.name = 'hdbscan_clusters')
      
      # hdbscan - UMAP
      print('UMAP reduction clusters...')
      hdbscan_res <- dbscan::hdbscan(integrated_data@reductions$umap@cell.embeddings,
                                     minPts = 50)
      integrated_data <- AddMetaData(integrated_data,
                                     metadata = hdbscan_res$cluster,
                                     col.name = 'hdbscan.UMAP_clusters')
      
      
    }else if(select_strategy == 'monocle3'){
      # monocle3 - PCA
      print('PCA reduction clusters...')
      cds <- SeuratWrappers::as.cell_data_set(integrated_data)
      cds <- monocle3::cluster_cells(cds, resolution=1e-3, reduction_method = 'PCA')
      
      integrated_data <- Seurat::AddMetaData(integrated_data,
                                             metadata = cds$ident,
                                             col.name = 'monocle3_clusters')
      # monocle3 - UMAP
      print('UMAP reduction clusters...')
      cds <- SeuratWrappers::as.cell_data_set(integrated_data)
      cds <- monocle3::cluster_cells(cds, resolution=1e-3, reduction_method = 'UMAP')
      
      integrated_data <- Seurat::AddMetaData(integrated_data,
                                             metadata = cds$ident,
                                             col.name = 'monocle3.UMAP_clusters')
      # release object
      rm(cds)
    }else if(select_strategy == 'scCCESS'){
      print('scCCESS employs its own dim reduction and clustering strategy...')
      res <- scCCESS::estimate_k(integrated_data@assays$integrated@data,
                                seed = 42, 
                                cluster_func = function(x,centers) { 
                                  set.seed(42);
                                  kmeans(x, centers)
                                },
                                criteria_method = "NMI",
                                krange = 5:15, ensemble_sizes = 10,
                                cores = parallel::detectCores()
      )
      
      res <- scCCESS::ensemble_cluster(integrated_data@assays$integrated@data,
                                       cluster_func = function(x) {
                                         set.seed(1)
                                         kmeans(x, centers = res$ngroups)
                                       }, 
                                       genes_as_rows = T,
                                       seed = 42,
                                       cores = parallel::detectCores(),
                                       ensemble_sizes = 10,
                                       scale = F, 
                                       batch_size = 64)
      
      
      integrated_data <- Seurat::AddMetaData(integrated_data,
                                             metadata = res,
                                             col.name = 'scCCESS_clusters')
    } else{
      print('Unknown strategy, skipping.')
    }
    
  }
  
  return(integrated_data)
}


#' Helper function for running FindMarkers on select ids.
#' @param integrated_data A seurat object
#' @param id The cluster id to find DEG in.
#' @return a seurat object
find_markers_on_id <- function(integrated_data, ids = c('seurat_clusters', 'seurat.UMAP_clusters', 'hdbscan_clusters', 'hdbscan.UMAP_clusters', 'monocle3_clusters', 'monocle3.UMAP_clusters','scCCESS_clusters')){

  DEG_markers <- list()
  for (id in c(ids)){
    print(paste('Running for cluster: ', id))
    integrated_data <- SetIdent(integrated_data, value = id) 
    DEG_markers[[id]] <- FindAllMarkers(integrated_data)
  }
  return(DEG_markers)
}

#' Density scoring based on 
#' @param integrated_data A seurat object
#' @param ct_markers A list of cell types and their corresponding markers
#' @param DEG_markers Markers identified by FindMarkers on said clustering, through find_markers_on_id
#' @param annotation_name Markers identified by FindMarkers on said clustering, through find_markers_on_id
#' @param clustering Clustering annotation to use for density annotation
#' @param cutoff A quality control cutoff to use for discerning cell types following scoring
#' @return a seurat object

density_score <- function(integrated_data, ct_markers, DEG_markers_set, annotation_name, clusterings = c("seurat_clusters", "seurat.UMAP_clusters", "hdbscan_clusters", "hdbscan.UMAP_clusters", "monocle3_clusters", "monocle3.UMAP_clusters", "scCCESS_clusters"), cutoffs = c('0.5', 'mean', 'none')){
  
  for (clustering in clusterings){
    
    if (!(clustering %in% colnames(integrated_data@meta.data))){
      print(paste('Clustering', clustering, 'is not computed, skipping...'))
    } else{
  
      integrated_data <- SetIdent(integrated_data, value = clustering) 
      
      DEG_markers <- DEG_markers_set[[clustering]]
      
      for (cutoff in cutoffs){
        # create a matrix to store the density scores in
        cell_type_satisfaction_matrix_clusters <- matrix(0, nrow = length(markers), ncol = length(levels(integrated_data)))
        rownames(cell_type_satisfaction_matrix_clusters) <- names(markers)
        
        cluster_i <- 0
        for (cluster_i in seq_along(levels(integrated_data)) - 1){
          
          curr_markers <- DEG_markers[DEG_markers$cluster == cluster_i, ] # get markers for current cluster
          
          curr_markers <- curr_markers[curr_markers$avg_log2FC > 1, ] # get only very up-regulated genes
          
          for (marker_set in names(markers)){
            
            curr_cell_type_set <- markers[[marker_set]] # get current markers
            
            top_genes <- curr_markers$gene # up-reguloted genes
            
            intersecting_genes <- intersect(curr_cell_type_set, top_genes) # present in the current marker
            
            cell_type_satisfaction_matrix_clusters[marker_set, cluster_i + 1] <- sum(curr_markers[which(curr_markers$gene %in% intersecting_genes), 'avg_log2FC']) / length(curr_cell_type_set)
            
            # prevent size of 1 cell types to have a high score by scaling their score by 0.8.
            if ((length(curr_cell_type_set) <= 1)){
              cell_type_satisfaction_matrix_clusters[marker_set, cluster_i + 1] <- cell_type_satisfaction_matrix_clusters[marker_set, cluster_i + 1] * 0.8
            }
          }
        }
        
        
        # get maximum scores for each cluster and find their identity
        max_intersection_values <- apply(cell_type_satisfaction_matrix_clusters, 2, max)
        
        top_cluster_res.manual <- c()
        for (i in seq_along(max_intersection_values)){
          res <- names(which(cell_type_satisfaction_matrix_clusters[, i] == max_intersection_values[i]))
          
          if (length(res) > 1){
            res <- 'Undecided'
          }
          if (cutoff != 'none'){
            if ((cutoff == 'mean') & (max_intersection_values[i] < mean(max_intersection_values))){
              res <- 'Undecided'
            }
            else if ((cutoff == '0.5') & (max_intersection_values[i] < 0.5)){
              res <- 'Undecided'
            }
          }
          top_cluster_res.manual <- c(top_cluster_res.manual, res)
        }
        
        names(top_cluster_res.manual) <- 1:length(top_cluster_res.manual) - 1
        new_metadata <- rep(0, ncol(integrated_data))
        
        # match each cluster with identified cell types
        for (i in seq_along(top_cluster_res.manual)){
          new_metadata[(integrated_data@meta.data[[clustering]] == names(top_cluster_res.manual[i]))] <- top_cluster_res.manual[i]
        }
        
        # for cell types with no markers
        new_metadata[new_metadata == '0'] <- 'Undecided'
        
        metadata_name <- paste0(annotation_name, '_', clustering, '_', cutoff)
        metadata_name <- gsub(" ", ".", metadata_name)
        metadata_name <- gsub("/", ".", metadata_name)
        
        # add annotation
        integrated_data <- AddMetaData(integrated_data,
                                       metadata = new_metadata,
                                       col.name = metadata_name)
      }
    }
  }
  
  
  
  return(integrated_data)
}

#' Helper function for loading results from python
#' @param integrated_data A seurat object
#' @param id The cluster id to find DEG in.
#' @return a seurat object
load_dnn_results <- function(integrated_data, new_annotations){
  
  # load up all the results and add them
  for (annotation in colnames(new_annotations)){
    if (!(new_annotation %in% colnames(integrated_data@meta.data))){
      integrated_data <- AddMetaData(integrated_data,
                                     metadata = new_DGCyTOF_annotations[,new_annotation],
                                     col.name = new_annotation)
    }
  }
  
  return(integrated_data)
}
 

# 
# 
# #################### Export To Python for DGCyTOF (Loom)  ###################
# 
# SaveLoom(integrated_data, filename = file.path(save_loc, "integrated_data.loom"))
# 
# #################### /Export To Python DGCyTOF (Loom)  ###################
# 
