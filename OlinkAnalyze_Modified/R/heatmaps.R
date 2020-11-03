#' Create heatmap for olink proteomics data
#' 
#' @return 
#' pheatmap object
#' 
#' @export
#' @import dplyr stringr tidyr pheatmap

olink_heatmap <- function(
  df, 
  GeneIDs = NULL, 
  SampleIDs = NULL, 
  verbose=TRUE,
  label_rows = FALSE, 
  label_cols = FALSE,
  row_clusters = NA,
  col_clusters = NA
) {
  
  # filter based on genes/samples
  if (!is.null(GeneIDs)) {
    df <- filter(df, GeneID %in% GeneIDs)
  }
  
  if (!is.null(SampleIDs)) {
    df <- filter(df, SampleID %in% SampleIDs)
  }
  
  #wide format
  df_wide <- df %>% 
    select(SampleID, Index, UniqueGeneID, NPX, case.control, WHO.severity) %>% 
    filter(!is.na(NPX)) %>% 
    spread(UniqueGeneID, NPX)
  
  percent_missingness <- colSums(is.na(df_wide[, -c(1:4)]))/nrow(df_wide)
  
  # assays with missingness > 20% are dropped from the PCA
  PERCENT_CUTOFF <- 0.2
  
  if(any(percent_missingness > PERCENT_CUTOFF)){
    
    removed_assays_index <- which(percent_missingness > PERCENT_CUTOFF)
    percent_missingness <- percent_missingness[-removed_assays_index]
    
    removed_assays_index <- removed_assays_index + 4
    removed_assays <- colnames(df_wide)[removed_assays_index]
    
    df_wide <- df_wide[, -removed_assays_index]
    
    if(verbose){
      warning(paste0("There are ",
                     paste0(length(removed_assays)), 
                     " assay(s) dropped due to high missingness (>",
                     round(PERCENT_CUTOFF*100),
                     "%)."))
    }
  }
  
  #<= PERCENT_CUTOFF assays imputed
  
  if(any(percent_missingness <= PERCENT_CUTOFF & percent_missingness > 0)){
    
    imputed_assays_index <- which(percent_missingness <= PERCENT_CUTOFF & percent_missingness > 0)
    percent_missingness <- percent_missingness[-imputed_assays_index]
    
    imputed_assays_index <- imputed_assays_index + 4
    imputed_assays <- colnames(df_wide)[imputed_assays_index]
    
    df_wide <- df_wide %>%
      mutate_at(tidyselect::all_of(imputed_assays), 
                ~ifelse(is.na(.x), median(.x, na.rm = TRUE), .x))
    
    if(verbose){
      warning(paste0("There are ",
                     paste0(length(imputed_assays)), 
                     " assay(s) that were imputed by their medians."))
    }
  }
  
  if(!all(colSums(is.na(df_wide[, -c(1:4)])) == 0)){
    stop('Missingness imputation failed.')
  }
  
  # convert data to the format pheatmap wants
  ccontrol <- data.frame(case.control=df_wide$case.control, WHO.severity=df_wide$WHO.severity)
  rownames(ccontrol) <- df_wide$SampleID
  
  ppanel <- data.frame(panel = gsub("^[^_]*_", "", colnames(df_wide)))
  rownames(ppanel) <- colnames(df_wide)
  
  df_wide_matrix <- df_wide %>% 
    select(-Index, -case.control, -WHO.severity) %>%
    column_to_rownames('SampleID') %>%
    as.matrix
  
  hmap <- pheatmap::pheatmap(
    df_wide_matrix,
    show_rownames = label_rows,
    show_colnames = label_cols,
    annotation_row = ccontrol,
    annotation_col = ppanel,
    cutree_cols = col_clusters,
    cutree_rows = row_clusters,
    scale="row"
  )
  
  return(hmap)
}

#' Create corrplot for olink proteomics data
#' 
#' @param t
#' (logical) whether to transpose prior to generation of correlation matrix -
#' by default will be sample vs sample
#' 
#' 
#' @return 
#' pheatmap object
#' 
#' @export
#' @import dplyr stringr tidyr pheatmap

olink_corrplot <- function(
  df, 
  t=FALSE,
  plot_labels=FALSE, 
  GeneIDs=NULL,
  SampleIDs=NULL,
  verbose=TRUE,
  clusters=NA,
  additional_annotation=NULL
) {
  
  # filter based on genes/samples
  if (!is.null(GeneIDs)) {
    df <- filter(df, GeneID %in% GeneIDs)
  }
  
  if (!is.null(SampleIDs)) {
    df <- filter(df, SampleID %in% SampleIDs)
  }
  
  #wide format
  if (t) {
    df_wide <- df %>% 
      select(SampleID, Index, UniqueGeneID, NPX) %>% 
      filter(!is.na(NPX)) %>% 
      spread(UniqueGeneID, NPX)
    
    add_index <- 2
  } else {
    df_wide <- df %>% 
      select(SampleID, Index, UniqueGeneID, NPX, case.control, WHO.severity) %>% 
      filter(!is.na(NPX)) %>% 
      spread(UniqueGeneID, NPX)
    
    add_index <- 4
  }
  
  percent_missingness <- colSums(is.na(df_wide[, -c(1:add_index)]))/nrow(df_wide)
  
  # assays with missingness > 20% are dropped from the PCA
  PERCENT_CUTOFF <- 0.2
  
  if(any(percent_missingness > PERCENT_CUTOFF)){
    
    removed_assays_index <- which(percent_missingness > PERCENT_CUTOFF)
    percent_missingness <- percent_missingness[-removed_assays_index]
    
    removed_assays_index <- removed_assays_index + add_index
    removed_assays <- colnames(df_wide)[removed_assays_index]
    
    df_wide <- df_wide[, -removed_assays_index]
    
    if(verbose){
      warning(paste0("There are ",
                     paste0(length(removed_assays)), 
                     " assay(s) dropped due to high missingness (>",
                     round(PERCENT_CUTOFF*100),
                     "%)."))
    }
  }
  
  #<= PERCENT_CUTOFF assays imputed
  
  if(any(percent_missingness <= PERCENT_CUTOFF & percent_missingness > 0)){
    
    imputed_assays_index <- which(percent_missingness <= PERCENT_CUTOFF & percent_missingness > 0)
    percent_missingness <- percent_missingness[-imputed_assays_index]
    
    imputed_assays_index <- imputed_assays_index + add_index
    imputed_assays <- colnames(df_wide)[imputed_assays_index]
    
    df_wide <- df_wide %>%
      mutate_at(tidyselect::all_of(imputed_assays), 
                ~ifelse(is.na(.x), median(.x, na.rm = TRUE), .x))
    
    if(verbose){
      warning(paste0("There are ",
                     paste0(length(imputed_assays)), 
                     " assay(s) that were imputed by their medians."))
    }
  }
  
  if(!all(colSums(is.na(df_wide[, -c(1:add_index)])) == 0)){
    stop('Missingness imputation failed.')
  }
  
  if (t) {
    annotations <- data.frame(panel = gsub("^[^_]*_", "", colnames(df_wide)))
    rownames(annotations) <- colnames(df_wide)
    
    df_wide_matrix <- df_wide %>% 
      select(-Index) %>%
      column_to_rownames('SampleID') %>%
      as.matrix
    
  } else {
    annotations <- data.frame(case.control=df_wide$case.control, WHO.severity=df_wide$WHO.severity)
    rownames(annotations) <- df_wide$SampleID
    annotations$WHO.severity[is.na(annotations$WHO.severity)] <- "NEGATIVE"
    
    df_wide_matrix <- df_wide %>% 
      select(-Index, -case.control, -WHO.severity) %>%
      column_to_rownames('SampleID') %>%
      as.matrix
    
    df_wide_matrix <- t(df_wide_matrix)
  }
  
  # if additional annotation is available, convert it to correct format
  if (!is.null(additional_annotation)) {
    annotations$SampleID <- rownames(annotations)
    additional_annotation <- data.frame(SampleID=names(additional_annotation),
                                        kmeans=as.character(additional_annotation))
    
    annotations <- left_join(annotations, additional_annotation, by=c("SampleID" = "SampleID"))
    sID <- annotations$SampleID
    annotations <- data.frame(WHO.severity=annotations$WHO.severity,
                              kmeans=annotations$kmeans)
    rownames(annotations) <- sID
  }
  
  cor_matrix <- cor(df_wide_matrix)
  
  hmap <- pheatmap::pheatmap(
    cor_matrix,
    show_colnames = plot_labels,
    show_rownames = FALSE,
    annotation_col = annotations,
    annotation_row = annotations,
    cutree_cols = clusters,
    cutree_rows = clusters
  )
  
  return(hmap)
}
