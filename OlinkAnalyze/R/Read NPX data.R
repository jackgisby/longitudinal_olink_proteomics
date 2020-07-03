#' Function to read NPX data into long format
#'
#' Imports an NPX file exported from NPX Manager. 
#' No alterations to the output NPX Manager format is allowed.
#'
#' @param filename Path to file NPX Manager output file.
#' @return A tibble in long format.
#' @keywords NPX
#' @export
#' @examples \donttest{read_NPX("~/NPX data.xlsx")}
#' @import dplyr stringr tidyr

read_NPX <- function(filename, sample_manifest=NULL, tab=1){
  
  if (!is.null(sample_manifest)) {
    manifest <- read.csv(sample_manifest, sep="\t")
    colnames(manifest) <- sub(colnames(manifest), pattern="Sample.ID", replacement = "SampleID")
  }
  
  NORM_FLAG <-  F
  
  # hard coded number of lines to skip as it differs between each tab haha :)
  if (tab == 1) {
    skip_mod <- 0
  } else if (tab == 2 | tab == 3) {
    skip_mod <- -2
  } else if (tab == 4) {
    skip_mod <- 1
  } else if (tab == 5) {
    skip_mod <- -1
  }
  
  meta_dat <-  readxl::read_excel(filename, skip = 12 + skip_mod, n_max = 4,col_names = F,.name_repair="minimal", sheet=tab)
  
  NR_DEVIATIONS <- sum(stringr::str_detect(meta_dat[1,], 'QC Deviation from median'))
  print(meta_dat)
  control_index <- (stringr::str_detect(meta_dat[2,], 'Det Ctrl') |
                      stringr::str_detect(meta_dat[2,], 'Inc Ctrl 2') |
                      stringr::str_detect(meta_dat[2,], 'Inc Ctrl 1') |
                      stringr::str_detect(meta_dat[2,], 'Ext Ctrl'))
  
  meta_index <- which(is.na(control_index))
  
  for (i in meta_index) {
    control_index[i] <- FALSE
  }
  
  meta_dat[2, control_index] <- meta_dat[2, control_index]
  meta_dat[2, meta_index] <- meta_dat[1, meta_index]
  meta_dat[3, control_index] <- meta_dat[2, control_index]
  meta_dat[3, meta_index] <- meta_dat[1, meta_index]

  meta_dat[4,] <- meta_dat[3,]
  meta_dat[4,1] <- 'SampleID'
  NR_CONTROLS <- sum(control_index)
  
  nr_col<-ncol(meta_dat)
  names(meta_dat)<-as.character(1:nr_col)
  
  meta_dat<-meta_dat %>%
    rename(Name = "1")
  
  dat <- readxl::read_excel(filename, skip = 17 + skip_mod, col_names = F,.name_repair="minimal", col_types = c('text'), sheet=tab)
  
  nr_col<-ncol(dat)
  names(dat)<-as.character(1:nr_col)
  
  dat<-dat %>%
    rename(Name = "1")
  
  missfreq<-dat %>% filter(stringr::str_detect(Name, "Missing Data freq."))
  LOD<-dat %>% filter(stringr::str_detect(Name, "LOD"))
  norm_method <- dat %>% filter(stringr::str_detect(Name, "Normalization"))
  
  if(nrow(norm_method) == 0){
    dat <- dat[c(-1*(nrow(dat) - 2):nrow(dat)),]
  }else{
    dat <- dat[c(-1*(nrow(dat) - 3):nrow(dat)),]
    NORM_FLAG <- T
  }
  
  meta_dat<-rbind(meta_dat,missfreq,LOD,norm_method)
  nr_panel <- 1
  SampleID<-dat$Name
  
  Index_nr<-c(1:length(SampleID))
  
  panel_data<-list() ##NPX values to be stored
  QC_list<-list()    ##QC data
  meta_data_list<-list() ## meta data
  panel_list<-list()  ## combination of panel data and QC
  assay_name_list<-list()
  panel_list_long<-list()
  deviations_list <- list()
  
  BASE_INDEX <- 90
  
  if(NR_CONTROLS > 0){
    BASE_INDEX <- BASE_INDEX + NR_CONTROLS/nr_panel
  }
  
  i = 1
    
  panel_data[[i]]<-dat[,(2+((i-1)*BASE_INDEX)):((BASE_INDEX+1)+((i-1)*BASE_INDEX))]
  
  if(NR_DEVIATIONS == 0){
    
    QC_list[[i]]<-dat[,c((2+((nr_panel)*BASE_INDEX)+(i-1)),
                         (2+((nr_panel)*BASE_INDEX)+(i-1))+nr_panel)]
    
    meta_data_list[[i]]<-meta_dat[,c((2+((i-1)*BASE_INDEX)):((BASE_INDEX+1)+((i-1)*BASE_INDEX)),
                                     (2+((nr_panel)*BASE_INDEX)+(i-1)),
                                     (2+((nr_panel)*BASE_INDEX)+(i-1))+nr_panel)]
    
    
  }else{
    
    QC_list[[i]]<-dat[,c((2+((nr_panel)*BASE_INDEX)+(i-1)),
                         (2+((nr_panel)*BASE_INDEX)+(i-1))+nr_panel, 
                         (2+((nr_panel)*BASE_INDEX)+(i-1))+2*nr_panel+(i-1), 
                         (2+((nr_panel)*BASE_INDEX)+(i-1))+2*nr_panel+(i-1)+1)]
    
    meta_data_list[[i]]<-meta_dat[,c((2+((i-1)*BASE_INDEX)):((BASE_INDEX+1)+((i-1)*BASE_INDEX)),
                                     (2+((nr_panel)*BASE_INDEX)+(i-1)),
                                     (2+((nr_panel)*BASE_INDEX)+(i-1))+nr_panel,
                                     (2+((nr_panel)*BASE_INDEX)+(i-1))+2*nr_panel, 
                                     (2+((nr_panel)*BASE_INDEX)+(i-1))+3*nr_panel)]
    
    meta_data_list[[i]][4,(BASE_INDEX+3)] <- "QC Deviation Inc Ctrl"
    meta_data_list[[i]][4,(BASE_INDEX+4)] <- "QC Deviation Det Ctrl"
    
    
  }
  
  meta_data_list[[i]][4,(BASE_INDEX+1)] <- meta_data_list[[i]][2,(BASE_INDEX+1)]
  meta_data_list[[i]][4,(BASE_INDEX+2)] <- meta_data_list[[i]][2,(BASE_INDEX+2)]
  
  
  panel_list[[i]]<-cbind(panel_data[[i]],QC_list[[i]])
  
  colnames(panel_list[[i]]) <- unlist(meta_data_list[[i]][3,])
  
  
  panel_list[[i]][,c(-(BASE_INDEX+1),-(BASE_INDEX+2))] <- lapply(panel_list[[i]][,c(-(BASE_INDEX+1),-(BASE_INDEX+2))],
                                                                 function(x) as.numeric(stringr::str_replace_all(x, c('#' = '', ',' = '.', 'No Data' = NA))))
  
  
  assay_name_list[[i]]<-tibble(ID=c(t(meta_data_list[[i]][4,])),
                               Name=c(t(meta_data_list[[i]][2,])),
                               UniProt = c(t(meta_data_list[[i]][3,])),
                               Panel=c(t(meta_data_list[[i]][1,])),
                               MissingFreq=c(t(meta_data_list[[i]][5,])),
                               LOD = as.numeric(c(t(meta_data_list[[i]][6,]))))
  
  if(NORM_FLAG == T){
    assay_name_list[[i]] <- bind_cols(assay_name_list[[i]], 
                                      Normalization = c(t(meta_data_list[[i]][7,])))
  }
  
  #return(meta_data_list[[i]])
  
  panel_list_long[[i]] <- panel_list[[i]] %>%
    mutate(SampleID = SampleID) %>%
    mutate(Index = Index_nr) %>%
    gather(Assay, NPX, -SampleID,-"QC Warning",-"Plate ID",-Index,-matches("QC Deviation Inc Ctrl"), -matches("QC Deviation Det Ctrl")) %>%
    left_join(assay_name_list[[i]], by = c('Assay' = 'ID')) %>%
    select(SampleID,Index,Assay, UniProt, Name,MissingFreq,Panel,"Plate ID","QC Warning",LOD,NPX,matches("Normalization"), matches("QC Deviation Inc Ctrl"), matches("QC Deviation Det Ctrl")) %>%
    rename(PlateID ="Plate ID") %>%
    rename(QC_Warning = "QC Warning") %>%
    rename(OlinkID = Assay, Assay = Name)
    
  if (!is.null(sample_manifest)) {
    panel_list_long[[i]] <- panel_list_long[[i]] %>%
      left_join(manifest, by = c("SampleID", "SampleID"))
  }
  
  return(panel_list_long[[i]])
}


read_multitab_NPX <- function(filename, sample_manifest=NULL, num_tabs=1) {
  for (i in 1:num_tabs) {
    if (i == 1) {
      long_npx <- read_NPX(filename=filename, sample_manifest=sample_manifest)
    } else {
      long_npx <- rbind(long_npx,
                        read_NPX(filename=filename, 
                                 sample_manifest=sample_manifest,
                                 tab=i))
    }
  }
  
  return(long_npx)
}
