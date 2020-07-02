# O-Link Proteomics Raw Data Processing Functions 

#--------------- import
library(openxlsx)

#--------------- functions

#' @title
#' Process raw o-link proteomics datasets
#' 
#' @description 
#' Takes an excel file produced by o-link and extracts the underlying data 
#' matrix, as well as additional feature and sample-based metadata. 
#' 
#' @return 
#' Returns a list as follows:
#' \itemize{
#'     \item prot_matrix - Numeric data matrix containing features as rows and
#'     samples as columns.
#'     \item prot_metadata - Dataframe containing feature-based metadata. 
#'     Includes protein IDs and LOD information.
#'     \item sample_metadata - Dataframe containing sample-based metadata. 
#'     Includes sample IDs and chip information.
#' }
#' 
#' @author 
#' Jack Gisby
#' 
#' @export

process_proteomics_xlsx <- function(path_xlsx, sheet=1, max_above_lod=1) {
    prot_matrix <- read.xlsx(path_xlsx, colNames=FALSE, sheet=sheet)
    
    lod_remove <- which(as.numeric(prot_matrix[(nrow(prot_matrix) - 1),]) > max_above_lod)
    
    if (length(lod_remove) > 0) {
        prot_matrix <- prot_matrix[,-lod_remove]
    }
    
    # select row-based/feature metadata (first two and last two rows)
    prot_metadata <- prot_matrix[,-c(ncol(prot_matrix) - 1, ncol(prot_matrix))]
    prot_metadata <- data.frame(
        uniprot_id = as.character(prot_metadata[1,-1]),
        unique_id = as.character(prot_metadata[2,-1]),
        miss_data_freq = as.numeric(prot_metadata[(nrow(prot_metadata) - 1),-1]),
        LOD_max = as.numeric(prot_metadata[nrow(prot_metadata),-1])
    )
    
    # cut row-based metadata
    prot_matrix <- prot_matrix[-c(1, 2, nrow(prot_matrix) - 1, nrow(prot_matrix)),]
    
    # extract col-based/sample metadata
    sample_metadata <- data.frame(
        patient_id = prot_matrix[,1],
        flagged = prot_matrix[,(ncol(prot_matrix) - 1)],
        chip_name = prot_matrix[,ncol(prot_matrix)]
    )
    
    # remove col-based metadata and convert to numeric
    prot_matrix <- prot_matrix[,-c(1, ncol(prot_matrix) - 1, ncol(prot_matrix))]
    prot_matrix <- t(matrix(as.numeric(unlist(prot_matrix)),nrow=nrow(prot_matrix)))
    
    colnames(prot_matrix) <- sample_metadata$patient_id
    rownames(prot_matrix) <- 1:nrow(prot_matrix)
    
    return(list(prot_matrix=prot_matrix, 
                prot_metadata=prot_metadata, 
                sample_metadata=sample_metadata))
}

#' @title
#' Remove NA values from expression matrix
#' 
#' @description 
#' Replaces NA values within a data matrix with the row (feature) median
#' 
#' @return 
#' The same matrix with imputed values
#' 
#' @author 
#' Jack Gisby
#' 
#' @export

na_median_replace <- function(prot_matrix) {
    which_na <- which(is.na(prot_matrix), arr.ind=TRUE)
    
    # replace na values with row median (median across all proteins of same type)
    if (length(which_na) > 0) {
        for (i in seq_along(length(which_na[,1]))) {
            prot_matrix[which_na[i,1],which_na[i,2]] <- median(prot_matrix[which_na[i,1]])
        }
    }
    
    return(prot_matrix)
}
