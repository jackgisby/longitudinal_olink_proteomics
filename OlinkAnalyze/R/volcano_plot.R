#' Generate linear models and make volcano plot
#' 
#' @export
#' @import dplyr stringr tidyr ggfortify ggrepel gridExtra

olink_volcano_plot <- function(
  long,
  variable,            
  random,                    
  covariates = NULL
) {
  plasma_pvals <- olink_lmer(long, "case.control", random = random, covariates = covariates)
  
  plasma_models <- olink_lmer(long, "case.control", random = random, covariates = covariates, return.models = TRUE)
  
  fc <- sapply(plasma_models, function(lmer_model) {return(lmer_model@beta[2])})
  
  fc_df <- data.frame(UniqueGeneID=names(plasma_models), fc=fc)
  
  plasma_pvals <- left_join(plasma_pvals, fc_df)
  
  plasma_pvals$logp <- -log10(plasma_pvals$Adjusted_pval)
  
  plt <- ggplot(plasma_pvals, aes(x=fc, y=logp, colour=Adjusted_pval < 0.05)) +
    geom_point(alpha=0.5) +
    ylab("-log10(adj_p)") +
    xlab("Log2 Fold Change") +
    geom_text_repel(data=subset(plasma_pvals,((fc >= 1 | fc <= -0.6) & Adjusted_pval < 0.05) | logp > 6),
                    aes(fc, logp, label = UniqueGeneID), size = 3, color="steelblue") + 
    geom_hline(yintercept=6, alpha=0.5) + 
    geom_vline(xintercept=-0.6, alpha=0.5) +
    geom_vline(xintercept=1, alpha=0.5)
  
  return(plt)
}