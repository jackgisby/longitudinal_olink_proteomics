#'Function which performs a linear mixed model per protein
#'
#'Fits a linear mixed effects model for every protein (by UniqueGeneID) in every panel, using lmerTest::lmer and stats::anova.
#'The function handles both factor and numerical variables and/or covariates. \cr\cr
#'Samples that have no variable information or missing factor levels are automatically removed from the analysis (specified in a messsage if verbose = T).
#'Character columns in the input dataframe are automatically converted to factors (specified in a message if verbose = T). 
#'Numerical variables are not converted to factors. 
#'If a numerical variable is to be used as a factor, this conversion needs to be done on the dataframe before the function call. \cr\cr
#'Crossed analysis, i.e. A*B formula notation, is inferred from the variable argument in the following cases: \cr
#'\itemize{
#' \item c('A','B')
#' \item c('A:B') 
#' \item c('A:B', 'B') or c('A:B', 'A')
#'}
#'Inference is specified in a message if verbose = T. \cr 
#'For covariates, crossed analyses need to be specified explicity, i.e. two main effects will not be expaned with a c('A','B') notation. Main effects present in the variable takes precedence. \cr 
#'The random variable only takes main effect(s). \cr
#'The formula notation of the final model is specified in a message if verbose = T. \cr\cr
#'Output p-values are adjusted by stats::p.adjust according to the Benjamini-Hochberg method (“fdr”). 
#'Adjusted p-values are logically evaluated towards adjusted p-value<0.05. 
#'
#' @param df NPX data frame in long format with at least protein name (Assay), UniqueGeneID, UniProt, 1-2 variables with at least 2 levels.
#' @param variable Single character value or character array. 
#' Variable(s) to test. If length > 1, the included variable names will be used in crossed analyses .
#' Also takes ':'/'*' notation. 
#' @param outcome Character. The dependent variable. Default: NPX.
#' @param random Single character value or character array.
#' @param covariates Single character value or character array. Default: NULL.
#' Covariates to include. Takes ':'/'*' notation. Crossed analysis will not be inferred from main effects.
#' @param return.covariates Boolean. Default: False. Returns results for the covariates. Note: Adjusted p-values will be NA for the covariates.
#' @param verbose Boolean. Deafult: True. If information about removed samples, factor conversion and final model formula is to be printed to the console. 
#'
#' @return A tibble containing the results of fitting the linear mixed effects model to every protein by UniqueGeneID, ordered by ascending p-value. 
#' @export
#' @examples
#' \donttest{
#' # Results in model NPX~Time*Treatment+(1|Subject)+(1|Site)
#' lmer_results <- olink_lmer(df = npx_df,
#' variable=c("Time", 'Treatment'),
#' random = c('Subject', 'Site'))
#' }
#' @import dplyr stringr tidyr lmerTest purrr

olink_lmer <- function(df,                        
                       variable,                  
                       outcome="NPX",            
                       random,                    
                       covariates = NULL,         
                       return.covariates=F,
                       return.models=FALSE,
                       verbose=T,
                       reorder=NULL
) {  
  
  if(missing(df) | missing(variable) | missing(random)){
    stop('The df and variable and random arguments need to be specified.')
  }
  
  withCallingHandlers({
    
    #Allow for :/* notation in covariates
    variable <- gsub("\\*",":",variable)
    if(!is.null(covariates)) covariates <- gsub("\\*",":",covariates)
    
    add.main.effects <- NULL
    if(any(grepl(":",covariates))){
      tmp <- unlist(strsplit(covariates,":"))
      add.main.effects <- c(add.main.effects,setdiff(tmp,covariates))
      covariates <- union(covariates,add.main.effects)
    }
    if(any(grepl(":",variable))){
      tmp <- unlist(strsplit(variable,":"))
      add.main.effects <- c(add.main.effects,setdiff(tmp,variable))
      variable <- union(variable,unlist(strsplit(variable,":")))
      variable <- variable[!grepl(":",variable)]
    }
    
    #If variable is in both variable and covariate, keep it in variable or will get removed from final table
    covariates <- setdiff(covariates,variable)
    add.main.effects <- setdiff(add.main.effects, variable)
    
    #Variables to be checked
    variable_testers <- intersect(c(variable,covariates,random), names(df))
    
    ##Remove rows where variables or covariate is NA (cant include in analysis anyway)
    removed.sampleids <- NULL
    for(i in variable_testers){
      removed.sampleids <- unique(c(removed.sampleids,df$SampleID[is.na(df[[i]])]))
      df <- df[!is.na(df[[i]]),]
    }
    
    #Not testing assays that have all NA:s
    all_nas <- df  %>%
      group_by(UniqueGeneID) %>%
      summarise(n = n(), n_na = sum(is.na(!!rlang::ensym(outcome)))) %>%
      ungroup() %>%
      filter(n == n_na) %>%
      pull(UniqueGeneID)
    
    
    if(length(all_nas) > 0) {
      
      warning(paste0('The assays ',
                     paste(all_nas, collapse = ', '),
                     ' have only NA:s. They will not be tested.'),
              call. = F)
      
    }
    
    ##Convert character vars to factor
    converted.vars <- NULL
    num.vars <- NULL
    for(i in variable_testers){
      if (i == variable & !is.null(reorder)) {
        df[[i]] <- ordered(df[[i]], levels=reorder)
        converted.vars <- c(converted.vars,i)
      } else if(is.character(df[[i]])){
        df[[i]] <- factor(df[[i]])
        converted.vars <- c(converted.vars,i)
      } else if(is.numeric(df[[i]])){
        num.vars <- c(num.vars,i)
      }
    }
    
    
    #Not testing assays that have all NA:s in one level
    #Every sample needs to have a unique level of the factor
    
    nas_in_var <- character(0)
    
    if(!is.null(covariates)){
      factors_in_df <- names(df)[sapply(df, is.factor)] 
      single_fixed_effects <- c(variable,
                                intersect(covariates,
                                          factors_in_df))
    }else{
      single_fixed_effects <- variable
    }
    
    
    for(effect in single_fixed_effects){
      
      current_nas <- df %>%
        filter(!(UniqueGeneID %in% all_nas)) %>%
        group_by(UniqueGeneID, !!rlang::ensym(effect)) %>%
        summarise(n = n(), n_na = sum(is.na(!!rlang::ensym(outcome)))) %>%
        ungroup() %>%
        filter(n == n_na) %>%
        distinct(UniqueGeneID) %>%
        pull(UniqueGeneID)
      
      if(length(current_nas) > 0) {
        
        nas_in_var <- c(nas_in_var, current_nas)
        
        warning(paste0('The assay(s) ',
                       current_nas,
                       ' has only NA:s in atleast one level of ',
                       effect,
                       '. It will not be tested.'),
                call. = F)
      }
      
      number_of_samples_w_more_than_one_level <- df %>% 
        group_by(SampleID, Index) %>% 
        summarise(n_levels = n_distinct(!!rlang::ensym(effect), na.rm = T)) %>% 
        ungroup() %>% 
        filter(n_levels > 1) %>% 
        nrow(.)
      
      if (number_of_samples_w_more_than_one_level > 0) {
        stop(paste0("There are ", 
                    number_of_samples_w_more_than_one_level, 
                    " samples that do not have a unique level for the effect ",
                    effect, 
                    ". Only one level per sample is allowed."))
      }
      
      
    }
    
    if(!is.null(covariates)){
      formula_string <- paste0(outcome, "~", 
                               paste(variable,collapse="*"),
                               "+", 
                               paste(covariates, sep = '', collapse = '+'),
                               "+",
                               paste(paste0("(1|",random,")"),collapse="+"))
    }else{
      
      formula_string <- paste0(outcome, "~", paste(variable,collapse="*"),
                               "+",
                               paste(paste0("(1|",random,")"),collapse="+"))
    }
    
    #Get factors
    fact.vars <- sapply(variable_testers, function(x) is.factor(df[[x]]))
    fact.vars <- names(fact.vars)[fact.vars]
    
    
    #Print verbose message
    if(verbose){
      if(!is.null(add.main.effects) & length(add.main.effects) > 0){
        message("Missing main effects added to the model formula: ",
                paste(add.main.effects,collapse=", "))
      }
      if(!is.null(removed.sampleids) & length(removed.sampleids) >0){
        message("Samples removed due to missing variable or covariate levels: ",
                paste(removed.sampleids,collapse=", "))
      }
      if(!is.null(converted.vars)){
        message(paste0("Variables and covariates converted from character to factors: ",
                       paste(converted.vars,collapse = ", ")))
      }
      if(!is.null(num.vars)){
        message(paste0("Variables and covariates treated as numeric: ",
                       paste(num.vars,collapse = ", ")))
      }
      message("Linear mixed effects model fit to each assay: ",formula_string)
    }
    
    if(!is.null(covariates) & any(grepl(":", covariates))){
      covariate_filter_string <- covariates[str_detect(covariates, ':')]
      covariate_filter_string <- sub("(.*)\\:(.*)$", "\\2:\\1", covariate_filter_string)
      covariate_filter_string <- c(covariates, covariate_filter_string)
    
    }else{
      covariate_filter_string <- covariates
    }
    
    if (return.models) {
      #make LMM
      lmer_model<- df %>%
        filter(!(UniqueGeneID %in% all_nas)) %>%
        filter(!(UniqueGeneID %in% nas_in_var)) %>%
        group_by(Assay, UniqueGeneID, UniProt, Panel) %>%
        group_map(~single_lmer(data=.x, formula_string = formula_string))
      
      mod_names <- df %>%
        filter(!(UniqueGeneID %in% all_nas)) %>%
        filter(!(UniqueGeneID %in% nas_in_var)) %>%
        group_by(Assay, UniqueGeneID, UniProt, Panel) %>%
        group_data() %>% 
        select("UniqueGeneID")
      
      names(lmer_model) <- as.vector(mod_names$UniqueGeneID)
      
      return(lmer_model)
      
    } else {
      ##make LMM
      lmer_model<-df %>%
        filter(!(UniqueGeneID %in% all_nas)) %>%
        filter(!(UniqueGeneID %in% nas_in_var)) %>%
        group_by(Assay, UniqueGeneID, UniProt, Panel) %>%
        group_modify(~tidy(anova(single_lmer(data=.x, formula_string = formula_string)))) %>%
        ungroup() %>%
        mutate(covariates = term %in% covariate_filter_string) %>% 
        group_by(covariates) %>% 
        mutate(Adjusted_pval=p.adjust(p.value,method="fdr")) %>%
        mutate(Threshold  = ifelse(Adjusted_pval<0.05,"Significant","Non-significant")) %>%
        mutate(Adjusted_pval = ifelse(covariates,NA,Adjusted_pval),
               Threshold = ifelse(covariates,NA,Threshold)) %>% 
        ungroup() %>% 
        select(-covariates) %>% 
        arrange(p.value)
      
      if(return.covariates){
        return(lmer_model)
      } else{
        return(lmer_model %>% filter(!term%in%covariate_filter_string))
      }
    }
    
  }, warning = function(w) {
    if (grepl(x = w, pattern = glob2rx("*not recognized or transformed: NumDF, DenDF*")) |
        grepl(x = w, pattern = glob2rx("*contains implicit NA, consider using*"))){
      invokeRestart("muffleWarning") 
    }
  })
  
}

single_lmer <- function(data, formula_string) {
  
  out.model <- tryCatch(lmerTest::lmer(as.formula(formula_string),
                                       data=data,
                                       REML=F,
                                       control = lmerControl(check.conv.singular = "ignore")),
                        warning = function(w){
                          return(
                            lmerTest::lmer(as.formula(formula_string),
                                           data=data,
                                           REML=F,
                                           control=lmerControl(optimizer = "Nelder_Mead",
                                                               check.conv.singular = "ignore"))
                          )
                          
                        }
  )
  
  
  if(class(out.model)=="lmerModLmerTest") {
    return(out.model)
    
  } else {
    stop("Convergence issue not caught by single_lmer")
  }
}

#' simple mixed differential expression analysis plot
#' @export
#' @import ggplot2 dplyr

simple_mixed_de <- function(
  long, 
  case_control_to_remove=NA, 
  variable="case.control", 
  covariates=NULL,
  random = c("Individual.ID"),
  logp_label=6,
  fc_label=c(-0.5, 1),
  plot_xlab="Log2 Fold Change (P-N)",
  plot_title="plasma case/control",
  return_pvals=FALSE
) {
  
  if (!is.na(case_control_to_remove)) {
    long <- filter(long, case.control != case_control_to_remove)
  }
  
  cctrl_cov_pvals <- olink_lmer(long, variable, random=random, covariates=covariates)
  cctrl__cov_models <- olink_lmer(long, variable, random=random, covariates=covariates, return.models = TRUE)
  
  fc <- sapply(cctrl__cov_models, function(lmer_model) {return(lmer_model@beta[2])})
  
  fc_df <- data.frame(UniqueGeneID=names(cctrl__cov_models), fc=fc)
  cctrl_cov_pvals <- left_join(cctrl_cov_pvals, fc_df)
  cctrl_cov_pvals$logp <- -log10(cctrl_cov_pvals$Adjusted_pval)
  
  if (return_pvals) {
    return(cctrl_cov_pvals)
  }
  
  ggplot(cctrl_cov_pvals, aes(x=fc, y=logp, colour=factor(ifelse(Adjusted_pval < 0.05, "<0.05", "NS"), levels=c("NS", "<0.05")))) +
    geom_point(alpha=0.5) +
    ylab("-log10(Adjusted Pvalue)") +
    xlab(plot_xlab) +
    labs(colour = "Adjusted Pvalue") +
    ggtitle(plot_title) +
    geom_text_repel(data=subset(cctrl_cov_pvals,
                                logp > logp_label | 
                                  ((fc > fc_label[2] | fc < fc_label[1]) & Adjusted_pval < 0.05)),
                    aes(fc, logp, label = UniqueGeneID), size = 3, color="steelblue")
}
