#' kmeans
#' 
#' @export

olink_kmeans <- function (
  df, 
  k=4,
  x_val = 1, 
  y_val = 2, 
  verbose = T,
  plot=TRUE
) {
  #wide format
  df_wide_matrix <- median_impute(df)
  
  matrix_dist <- dist(df_wide_matrix)
  
  km <- kmeans(matrix_dist, k)
  
  df %>%
    left_join(data.frame(SampleID=names(km$cluster), km_clusts=as.character(km$cluster)), by=c("SampleID" = "SampleID")) %>%
    olink_pca_plot(color_g="km_clusts", y_val=y_val, x_val=x_val) %>%
    print()
  
  return(km)
}

#' consensus clustering
#' 
#' identify what k should be set to
#' 
#' @export

olink_consensus_clustering <- function(
  df,
  clusteralg="km",
  iters=25,
  reps=100
) {
  df_wide_matrix <- median_impute(df)
  
  return(M3C(t(df_wide_matrix), iters=iters, repsref=reps, repsreal=reps, clusteralg=clusteralg))
}

#' median imputation and transformation to wide matrix

median_impute <- function(
  df
) {
  verbose = TRUE
  
  df_wide <- df %>% 
    select(SampleID, Index, UniqueGeneID, NPX) %>% 
    filter(!is.na(NPX)) %>% 
    spread(UniqueGeneID, NPX)
  
  percent_missingness <- colSums(is.na(df_wide[, -c(1:2)]))/nrow(df_wide)
  
  # assays with missingness > ?% are dropped from the PCA
  PERCENT_CUTOFF <- 0.5
  
  if(any(percent_missingness > PERCENT_CUTOFF)){
    
    removed_assays_index <- which(percent_missingness > PERCENT_CUTOFF)
    percent_missingness <- percent_missingness[-removed_assays_index]
    
    removed_assays_index <- removed_assays_index + 2
    removed_assays <- colnames(df_wide)[removed_assays_index]
    
    df_wide <- df_wide[, -removed_assays_index]
    
    if(verbose){
      warning(paste0("There are ",
                     paste0(length(removed_assays)), 
                     " assay(s) dropped due to high missingness (>",
                     round(PERCENT_CUTOFF*100),
                     "%)."))
    }
    
    if(!is.null(loadings_list)){
      
      dropped_loadings <- intersect(removed_assays, 
                                    loadings_list)
      
      
      if(length(dropped_loadings) > 0){
        
        if(verbose){
          warning(paste0("The loading(s) ",
                         paste0(dropped_loadings, collapse=", "),
                         " from the loadings_list are dropped due to high missingness. "))
        }
        
        loadings_list <- setdiff(loadings_list, dropped_loadings)
        
        if(length(loadings_list) == 0){
          
          loadings_list <- NULL
          
        }
      }
      
    }
    
  }
  
  #<= PERCENT_CUTOFF assays imputed
  
  if(any(percent_missingness <= PERCENT_CUTOFF & percent_missingness > 0)){
    
    imputed_assays_index <- which(percent_missingness <= PERCENT_CUTOFF & percent_missingness > 0)
    percent_missingness <- percent_missingness[-imputed_assays_index]
    
    imputed_assays_index <- imputed_assays_index + 2
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
  
  if(!all(colSums(is.na(df_wide[, -c(1:2)])) == 0)){
    stop('Missingness imputation failed.')
  }
  
  df_wide %>% 
    select(-Index) %>%
    column_to_rownames('SampleID') %>%
    as.matrix() %>%
    return()
}
