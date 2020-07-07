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
#' @import dplyr stringr tidyr biomaRt

read_NPX <- function(filename, sample_manifest=NULL, pheno=NULL, tab=1){
  
  if (!is.null(sample_manifest)) {
    manifest <- read.csv(sample_manifest, sep="\t")
    colnames(manifest) <- sub(colnames(manifest), pattern="Sample.ID", replacement = "SampleID")
    
    if (!is.null(pheno)) {
      pheno_df <- read.csv(pheno, sep="\t")
      colnames(pheno_df) <- sub(colnames(pheno_df), pattern="Research.ID", replacement = "Individual.ID")
      
      manifest <- manifest %>%
        left_join(pheno_df)
    }
  }
  
  NORM_FLAG <-  F
  
  # hard coded number of lines to skip as it differs between each tab :)
  if (tab == 1) {
    skip_mod <- 0
    panel <- "CardMet"
    
  } else if (tab == 2) {
    skip_mod <- -2
    panel <- "CVD2"
    
  } else if (tab == 3) {
    skip_mod <- -2
    panel <- "CVD3"
    
  } else if (tab == 4) {
    skip_mod <- 1
    panel <- "Inf"
    
  } else if (tab == 5) {
    skip_mod <- -1
    panel <- "ImmResp"
    
  }
  
  meta_dat <-  readxl::read_excel(filename, skip = 12 + skip_mod, n_max = 4,col_names = F,.name_repair="minimal", sheet=tab)
  
  NR_DEVIATIONS <- sum(stringr::str_detect(meta_dat[1,], 'QC Deviation from median'))
  
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
  dat <- dat[-which(is.na(dat[,1]))[1],]

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
  
  gene_ids <- uniprot_to_gene(data.frame(UniProt=as.character(meta_dat[2,-1]), target=as.character(meta_dat[1,-1])))
  gene_ids <- data.frame(cbind(colnames(gene_ids), t(gene_ids)))
  names(gene_ids) <- names(meta_dat)
  
  meta_dat<-rbind(meta_dat,missfreq,LOD,norm_method,gene_ids)
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
                               UniProt = c(t(meta_data_list[[i]][2,])),
                               Panel=c(t(meta_data_list[[i]][1,])),
                               MissingFreq=c(t(meta_data_list[[i]][5,])),
                               LOD = as.numeric(c(t(meta_data_list[[i]][6,]))),
                               GeneID=c(t(meta_data_list[[i]][18,])))
  
  if(NORM_FLAG == T){
    assay_name_list[[i]] <- bind_cols(assay_name_list[[i]], 
                                      Normalization = c(t(meta_data_list[[i]][7,])))
  }

  panel_list_long[[i]] <- panel_list[[i]] %>%
    mutate(SampleID = SampleID) %>%
    mutate(Index = Index_nr) %>%
    gather(Assay, NPX, -SampleID,-"QC Warning",-"Plate ID",-Index,-matches("QC Deviation Inc Ctrl"), -matches("QC Deviation Det Ctrl")) %>%
    left_join(assay_name_list[[i]], by = c('Assay' = 'ID')) %>%
    dplyr::select(SampleID, Index, Assay, UniProt, GeneID, Name, MissingFreq, Panel, "Plate ID", "QC Warning", LOD, NPX, matches("Normalization"), matches("QC Deviation Inc Ctrl"), matches("QC Deviation Det Ctrl")) %>%
    rename(PlateID ="Plate ID") %>%
    rename(QC_Warning = "QC Warning") %>%
    rename(OlinkID = Assay, Assay = Name)
  
  if (!is.null(sample_manifest)) {
    panel_list_long[[i]] <- panel_list_long[[i]] %>%
      left_join(manifest)
  }
  
  panel_list_long[[i]]$panel_name <- panel
  
  return(panel_list_long[[i]])
}

#' Multitab NPX data parsing
#' 
#' @export
#' @import dplyr stringr tidyr

read_multitab_NPX <- function(filename, sample_manifest=NULL, pheno=NULL, num_tabs=1) {
  for (i in 1:num_tabs) {
    if (i == 1) {
      long_npx <- read_NPX(filename=filename, sample_manifest=sample_manifest, pheno=pheno)
    } else {
      long_npx <- rbind(long_npx,
                        read_NPX(filename=filename, 
                                 sample_manifest=sample_manifest,
                                 pheno=pheno,
                                 tab=i))
    }
  }
  
  return(long_npx)
}

#' Map protein IDs to gene IDs. 
#' 
#' @param df meta_dat dataframe
#' @param tab (optional) tab number
#' 
#' @return Gene IDs

uniprot_to_gene <- function(df) {
  df$uniprot <- df$UniProt
  
  # data input error by Olink 'o' instead of 'O'
  df$target <- gsub("IL-2oRA", "IL-20RA", df$target)
  
  # turns out TWEAK labelled with an out of date or inferior UP id
  if (TRUE %in% grepl('TWEAK', df$target)) {
    df$uniprot[grep('TWEAK', df$target)] <- "O43508"
  }
  
  # no uniprot id for NT-proBNP - there is for proBNP
  if ("NT-proBNP" %in% df$target) {
    df[df$target == "NT-proBNP",]$uniprot <- "P16860"
  }
  # Clean up 1. identify bad entries: Olink have made some errors and Excel import causes some problems
  
  # clean whitespace
  df$uniprot <- gsub('\\\r\\\n', ";", df$uniprot, ignore.case = F)
  df$uniprot <- gsub(', |,', ";", df$uniprot, ignore.case = F)
  df$uniprot <- gsub("[[:space:]]", "", df$uniprot)
  
  # Clean up 2. '-' represents isoform notation eg O43521-2
  
  df$uniprot.isoform <- NA
  
  df$uniprot.isoform[grep('-', df$uniprot)] <- grep('-', df$uniprot, v=T)
  df$uniprot <- gsub("-[0-9]$", "", df$uniprot)

  # Special circumstances 2:two ids for protein complex eg IL12A-IL12B
  # uniprot ids sep by ';'
  # df[grep(";", df$uniprot), ]
  
  df$multiple.proteins <- FALSE
  df$multiple.proteins[grep(";", df$uniprot)] <- TRUE
  
  df$protein.1 <- df$uniprot  
  df$protein.2 <- NA
  
  df$protein.2[which(df$multiple.proteins==T)] <- str_extract(string=df$uniprot, pattern = ";[A-Z0-9]+")[which(df$multiple.proteins==T)]
  df$protein.2 <- gsub("^;", "", df$protein.2)
  
  df$protein.1[which(df$multiple.proteins==T)] <- str_extract(string=df$uniprot, pattern = "^[A-Z0-9]+;")[which(df$multiple.proteins==T)]
  df$protein.1 <- gsub(";$", "", df$protein.1)
  
  # where there are 2 uniprot ids (eg protein complex) the uniprot ids are not always in consistent order
  # lets make them in consistent alphabetical order
  df$uniprot.ordered <- NA
  df$uniprot.ordered[!df$multiple.proteins] <- df$uniprot[!df$multiple.proteins]
  
  alphabetize.up <- function(x) {
    if( !inherits(x, what='data.frame')){
      stop('argument must be a data.frame')
    }
    y <- paste(sort(c( x$protein.1, x$protein.2)),collapse=';')
    y
  }
  
  inds <- which(df$multiple.proteins)
  
  for (i in inds){
    df$uniprot.ordered[i] <- alphabetize.up(df[i,]) 
  }
  
  #annoying that p1 and p2 are arbitrary: now we've ordered things let's start over on this front
  df$protein.1 <- df$protein.2 <-NULL
  
  # now repeat the exercise for p1 and p2 using the alphabetized concatenation
  
  df$protein.1 <- df$uniprot.ordered  
  df$protein.2 <- NA
  
  df$protein.2[which(df$multiple.proteins==T)] <- str_extract(string=df$uniprot.ordered, pattern = ";[A-Z0-9]+")[which(df$multiple.proteins==T)]
  df$protein.2 <- gsub("^;", "", df$protein.2)
  
  df$protein.1[which(df$multiple.proteins==T)] <- str_extract(string=df$uniprot.ordered, pattern = "^[A-Z0-9]+;")[which(df$multiple.proteins==T)]
  df$protein.1 <- gsub(";$", "", df$protein.1)
  
  # col to identify dup proteins and which panels
  
  dup.prots <- union( which( duplicated(df$uniprot.ordered)), which( duplicated(df$uniprot.ordered, fromLast = T)) )
  df$prot.on.multiple.panel <- FALSE
  df$prot.on.multiple.panel[dup.prots] <- TRUE
  
  df$panels.with.prot <- NA
  
  
  tmp.list <- split( df[dup.prots,], f=df$uniprot.ordered[dup.prots] )
  
  mylist <- lapply(tmp.list, FUN = function(x) paste( as.character(x$panel), collapse=";" ) )
  
  for (i in dup.prots){
    uprot <- df$uniprot.ordered[i]
    df[i, "panels.with.prot"] <- mylist[[uprot]]
  }
  
  #--------------------- Gene symbol annotation ---------------------#
  
  # matching to gene symbols: do for p1 and p2
  
  #ensembl <- biomaRt::useMart(biomart="ensembl",
  #                   dataset="hsapiens_gene_ensembl",
  #                   host='http://jul2018.archive.ensembl.org')

  #filters <- biomaRt::listFilters(ensembl)

  #x <- biomaRt::getBM(attributes = c('uniprotswissprot', 'hgnc_symbol', 'entrezgene', 'chromosome_name'),
  #             filters = 'uniprotswissprot',
  #             values = df$protein.1,
  #             mart = ensembl)

  # some UP ids not found by BioMart: turns out we have outdated IDs
  #df[which(!df$protein.1 %in% x$uniprotswissprot),]
  
  #--------------------- Try an archived version of Ensembl ---------------------#
  
  # find urls for old ensembl versions
  biomaRt::listEnsemblArchives()
  
  # hg19/GRCh37
  ensembl.hg19 <- biomaRt::useMart(biomart= "ENSEMBL_MART_ENSEMBL",
                          dataset="hsapiens_gene_ensembl",
                          host = 'http://grch37.ensembl.org')
  
  # note attribute names differ in the older release
  gene.pos <- biomaRt::getBM(attributes = c('uniprotswissprot', 'hgnc_symbol', # 'entrezgene',
                                   'chromosome_name', 'start_position', 'end_position'), 
                    filters = 'uniprotswissprot', 
                    values = unique(df$protein.1), 
                    mart = ensembl.hg19)
  
  # P0DOY2 / IGLC2 is present within db but swissprot id is NA
  missing_to_add <- biomaRt::getBM(
    attributes = c('uniprotswissprot', 'hgnc_symbol', # 'entrezgene', 
                   'chromosome_name', 'start_position', 'end_position'), 
    filters = 'hgnc_symbol', 
    values = "IGLC2", 
    mart = ensembl.hg19)
  
  missing_to_add[missing_to_add$hgnc_symbol == "IGLC2",]$uniprotswissprot <- "P0DOY2"
  
  # same for hgnc_symbol
  missing_to_add <- rbind(missing_to_add, data.frame(biomaRt::getBM(
    attributes = c('uniprotswissprot', 'hgnc_symbol', # 'entrezgene', 
                   'chromosome_name', 'start_position', 'end_position'), 
    filters = 'hgnc_symbol', 
    values = "PECAM1", 
    mart = ensembl.hg19)))
  
  missing_to_add[missing_to_add$hgnc_symbol == "PECAM1",]$uniprotswissprot <- "P16284"
  
  # same for CDSN
  missing_to_add <- rbind(missing_to_add, data.frame(biomaRt::getBM(
    attributes = c('uniprotswissprot', 'hgnc_symbol', # 'entrezgene', 
                   'chromosome_name', 'start_position', 'end_position'), 
    filters = 'hgnc_symbol', 
    values = "CDSN", 
    mart = ensembl.hg19)))
  
  missing_to_add[missing_to_add$hgnc_symbol == "CDSN",]$uniprotswissprot <- "Q15517"

  gene.pos <- rbind(gene.pos, missing_to_add)

  # there are some duplicated genes
  dup.ind <- union( which(duplicated(gene.pos$hgnc_symbol)),
                    which(duplicated(gene.pos$hgnc_symbol, fromLast = T))
  )
  
  # strange chr names
  
  strange.ind <- which(!gene.pos$chromosome_name %in% c(1:22, 'X', 'Y'))
  
  to.cut <- intersect(dup.ind, strange.ind)
  
  gene.pos2 <- gene.pos[-to.cut,]

  #-------------------------------------------------------------------------------#
  # some proteins map to multiple genes
  map_ids <- vector("integer", nrow(df)) 
  for (i in 1:nrow(df)) {
    map_ids[i] <- which(gene.pos2$uniprotswissprot == df$protein.1[i])[1]
  }
  gene.pos2 <- gene.pos2[map_ids,]
  
  df2 <- left_join(df, gene.pos2, by=c("protein.1" = "uniprotswissprot"))

  return(df2)
}
