library(Seurat)
library(scales)
library(ggplot2)
library(devtools)
library(data.table)
library(dplyr)
library(future)
library(glue)
library(tibble)
library(pheatmap)
library(grid)
library(shiny)
library(shinythemes)
library(shinyWidgets)
library(shinyvalidate)
library(shinyalert)
library(shinydashboard)
library(shinyjs)
library(DT)
library(shinydashboardPlus)
library(ggplot2)
library(markdown)
library(ggthemes)
library(SeuratObject)
library(SeuratData)
library(clustree)


# Module containing functions for sc-pipeline app 

# Input/ config section 
load_seurat_obj <- function(path){
  errors <- c()
  
  #extension check 
  if(!tolower(tools::file_ext(path)) == "rds"){
    errors <- c(errors, "Invalid rds file.")
    return(errors)
  }
  
  tryCatch(
    {
      seurat <- readRDS(path)
    },
    error = function(e){
      errors <- c(errors, "Invalid rds file.")
      return(errors)
    }
  )
  
  if (!inherits(seurat, "Seurat")){
    errors <- c(errors, "File is not Seurat object.")
  }
  return(seurat)
}


# Pre-processing section. 

  ## Define variables based on input ids 
  preprocess <- function(obj, mt,nfU,nfL,assay,sct_features,pcs){
    withProgress(message = "Preprocessing data", value = 0,{
    
    ## Filtering subsection
    
    ### Calculate and filter %MT
    incProgress(0.25,message = "Filtering MT")    
    #message(class(obj)) <- for troubleshooting! - DM
    nfU <- as.integer(nfU)
    nfL <- as.integer(nfL)
    mt <- as.integer(mt)
    pcs <- as.integer(pcs)
    
    obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
    obj <- subset(obj, subset = percent.mt < mt)
    incProgress(0.25,message = "Filtering Features")
    
    ### nfeature filter
    message("Filtering nFeatures")
    x <- print(paste0("nFeature_", assay))
    obj <- subset(obj, nFeature_RNA < nfU & nFeature_RNA > nfL)
    
    incProgress(0.25,message = "Normalizing Data: This may take a while")
    
  ## Normalization, Scaling, SCT Assay, Variable Features (NOTE: eventually allow other selection methods - DM)
    message("Normalizing Data: This can take a while.")
    obj <- NormalizeData(obj)
    obj <- ScaleData(obj,verbose = F)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = sct_features)
    obj <- SCTransform(obj, assay = assay, verbose = F)
    obj <- RunPCA(obj, assay = assay, verbose = F)
    
    incProgress(0.25,message = "Clustering Data")
    
  ## Preliminary clustering 
    message("Clustering")
    obj <- FindNeighbors(obj, reduction = "pca", dims = 1:pcs, verbose = F)
    #message(names(obj)) <- for troubleshooting! - DM 
  
    obj <- FindClusters(obj, verbose = F, graph.name = "RNA_snn")
    obj <- RunUMAP(obj, reduction = 'pca', dims = 1:pcs, verbose = F)
    processed_obj <- obj
    message("Done!")
    return(processed_obj)
    })
  }
  
  preprocess_fields <- function(mt,nfU,nfL,assay,sct_features,pcs){
    
    blankfields <- c()
    if(mt == ""){
      blankfields <- c(blankfields, "Mitochondrial Filter")
    }
    if(nfL == ""){
      blankfields <- c(blankfields, "Minimum features per cell")
    }
    if(nfU == ""){
      blankfields <- c(blankfields, "Maximum features per cell")
    }
    if(sct_features == ""){
      blankfields <- c(blankfields, "SCT Features")
    }
    if(pcs == ""){
      blankfields <- c(blankfields, "Principal Components")
    }
    
    return(blankfields)
  }
  
# Integration Section (Currently supports merging 2 objects)
 integrate_reference <- function(source,reference){
     withProgress(message = "Integrating Reference Data", value = 0,{
     DefaultAssay(source) <- "RNA" 
     
     incProgress(0.33,message = "Defining Anchorset: This may take a while.")
     object.list <- c(reference,source)
     object.anchors <<- FindIntegrationAnchors(object.list = object.list, dims = 1:30)
     
     incProgress(0.33,message = "Integrating: This may take a while.")
     object.integrated <- IntegrateData(anchorset = object.anchors, dims = 1:30)
     
     incProgress(0.165, message = "Processing Integrated Assay")
     DefaultAssay(object.integrated) <- "integrated"
     object.integrated <- ScaleData(object.integrated, verbose = FALSE)
     object.integrated <- RunPCA(object.integrated, npcs = 30, verbose = FALSE)
     incProgress(0.165, message = "Clustering Integrated Object")
     object.integrated <- FindNeighbors(object.integrated, reduction = "pca", dims = 1:pcs, verbose = F)
     object.integrated <- FindClusters(object.integrated, verbose = F, graph.name = "integrated_snn")
     object.integrated <- RunUMAP(object.integrated, reduction = "pca", dims = 1:30, verbose = FALSE)
     message("Done!")
     return(object.integrated)
     })
   }
 
 
 # Annotation Section
 
 ## Majority Voting: Inspired by CellTypist (not America)
 democratic_annotation <- function(object){
   for(cluster in unique(object$seurat_clusters)){
     print(cluster) 
     cluster_cells <- object@meta.data %>% filter(seurat_clusters == cluster)
     cluster_annotations <- table(cluster_cells$seurat_annotations)
     print(which.max(cluster_annotations))
   }
 }
 
 # Differential Expression Section 
 
 marker_analysis <- function(object){
   markers <- FindAllMarkers(object, min.pct = 0.25,logfc.threshold = 0.25)
   markers %>%
     group_by(cluster) %>%
     slice_max(n = 10, order_by = avg_log2FC)
   
 }

 # New clustree based resolution optimizer
clustree_markers <- function(){
  
}
 
