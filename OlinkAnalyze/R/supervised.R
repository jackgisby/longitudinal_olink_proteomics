#' run random forests
#' @export
#' @import caret randomForestExplainer

run_rf <- function(long, variable="grouped_severity", sampling=NULL,
                   selected_features=NULL) {
    
    # convert data from long format to a matrix
    prot_matrix <- long %>% 
        dplyr::select(SampleID, GeneID, NPX, variable) %>% 
        filter(!is.na(NPX)) %>% 
        spread(GeneID, NPX) %>%
        column_to_rownames('SampleID')
    
    # select specific features for model
    if (!(is.null(selected_features))) {
        prot_matrix <- prot_matrix[,colnames(prot_matrix) %in% selected_features
                                   | colnames(prot_matrix) == variable]
    }
    
    # train model using caret
    rf <- train(as.formula(paste0(variable, " ~ .")), 
                data = prot_matrix, 
                localImp = TRUE,
                method = "rf", 
                metric = "Kappa",
                proximity=TRUE,
                preProcess = c("knnImpute", "scale", "center"),
                na.action = na.pass,
                tuneGrid = data.frame(mtry = floor(sqrt(ncol(prot_matrix) - 1))),
                trControl = trainControl(method = "repeatedcv", p = 0.7, number=10, repeats=10, sampling=sampling))
    
    rf$finalModel$call$formula <- eval(rf$call$form)
    print(rf)
    print(rf$finalModel)
    
    # plot feature importance (ntrees, minimal depth)
    print(plot_min_depth_distribution(min_depth_distribution(rf$finalModel), k = 15) +
              theme_minimal() +
              theme(title=element_blank(), panel.grid.major = element_blank()))
    
    return(rf)
}

#' plot random forests
#' @export
#' @import caret randomForestExplainer
#' 
plot_rf <- function(rf,
                    no_of_trees_label=40, accuracy_decrease_label=0.004, return_imp_frame=FALSE) {
    
    # get feature importances
    imp_frame <- measure_importance(rf$finalModel)
    
    # adjust p values
    imp_frame$p_adj <- p.adjust(imp_frame$p_value, "fdr")
    imp_frame$p_s <- factor(ifelse(imp_frame$p_adj < 0.05, "<0.05", "NS"), levels=c("NS", "<0.05"))
    
    if (return_imp_frame) {
        return(imp_frame)
    }
    
    # secondary feature importance plot
    imp_plot <- ggplot(imp_frame, aes(no_of_trees, accuracy_decrease, col=p_s)) + 
        theme_pubr() +
        xlab("Number of Trees") + ylab("Accuracy Decrease") +
        geom_point() + 
        labs(col = "PValue") +
        geom_text_repel(data=subset(imp_frame, no_of_trees > no_of_trees_label | accuracy_decrease_label > 0.004),
                        aes(no_of_trees, accuracy_decrease, label = variable), size = 3, color="steelblue")
    
    print(imp_plot)
    # print(plot_importance_ggpairs(imp_frame, c("no_of_trees", "accuracy_decrease", "p_adj", "mean_min_depth", "gini_decrease")) + 
    #           theme_pubr(border = TRUE))
    
    return(important_variables(imp_frame, k = 30, measures = c("mean_min_depth", "no_of_trees")))
}

#' plot random forest interactions
#' @export
#' @import caret randomForestExplainer

plot_rf_interactions <- function(long, rf, vars, variable="grouped_severity", 
                                 num_interactions=15, make_the_plot=TRUE) {
    
    # convert from long format to matrix
    prot_matrix <- long %>% 
        dplyr::select(SampleID, GeneID, NPX, variable) %>% 
        filter(!is.na(NPX)) %>% 
        spread(GeneID, NPX) %>%
        column_to_rownames('SampleID')
    
    if (make_the_plot) {
        # lots of warnings produced by function
        suppressWarnings(interactions_frame <- min_depth_interactions(rf$finalModel, vars))
        
        # look at commonly "interacting" featuress
        inter_plot <- plot_min_depth_interactions(interactions_frame, k=num_interactions) + 
            ylab("Mean minimum depth")
        
        all_ids_to_remove <- vector("logical", length(inter_plot$data$interaction))
        for (repeat_inter in unique(inter_plot$data$interaction[duplicated(inter_plot$data$interaction)])) {
            repeat_ids <- which(inter_plot$data$interaction == repeat_inter)
            selected_id <- which.min(inter_plot$data$mean_min_depth[repeat_ids])
            all_ids_to_remove[repeat_ids[-selected_id]] <- TRUE
        }
        
        if (any(all_ids_to_remove)) {
            inter_plot$data <- inter_plot$data[-which(all_ids_to_remove),]
        }
        
        inter_plot$layers <- NULL
        
        inter_plot <- inter_plot +
            geom_bar(stat="identity", aes(reorder(interaction, -occurrences), mean_min_depth)) +
            # geom_pointrange(stat = "identity",
            #                 aes(y=uncond_mean_min_depth, 
            #                     ymin = pmin(mean_min_depth, uncond_mean_min_depth), 
            #                     ymax = pmax(mean_min_depth, uncond_mean_min_depth))) +
            scale_fill_viridis(aes(fill=occurrences)) +
            theme_minimal() +
            xlab("Interaction") +
            guides(fill=guide_colourbar(title="Occurrences")) +
            theme(panel.grid.major.y = element_blank()) +
            coord_flip()
        
        print(inter_plot)
    }
    
    # generate imputed data matrix
    x <- prot_matrix %>%
        preProcess(method = c("knnImpute", "scale", "center")) %>%
        predict(prot_matrix %>% select(-variable))
    
    x[[variable]] <- prot_matrix[[variable]]
    
    return(x)
}

#' plot the interaction between two variables according to rf predictions
#' @export
#' @import randomForestExplainer reshape2
plot_interactions <- function(final_rf, x, var1, var2, grid=300, virtis_option="D") {
    
    inter_plot <- plot_predict_interaction(final_rf, 
                                           x, 
                                           var1, var2, 
                                           grid=grid)
    
    var1_min <- min(inter_plot$data[[var1]])
    var1_max <- max(inter_plot$data[[var1]])
    var1_breaks <- seq(var1_min, var1_max, (var1_max - var1_min) / 30)
    var2_min <- min(inter_plot$data[[var2]])
    var2_max <- max(inter_plot$data[[var2]])
    var2_breaks <- seq(var2_min, var2_max, (var2_max - var2_min) / 30)
    
    a <- melt(tapply(inter_plot$data$prediction,
                     list(x=cut(inter_plot$data[[var1]], labels=FALSE, breaks=var1_breaks), 
                          y=cut(inter_plot$data[[var2]], labels=FALSE, breaks=var2_breaks)),
                     mean))
    
    colnames(a) <- c("var1", "var2", "Prediction")
    
    a$var1 <- sapply(a$var1, function(x) {
        return(mean(var1_breaks[x], var1_breaks[x+1]))
    })
    
    a$var2 <- sapply(a$var2, function(x) {
        return(mean(var2_breaks[x], var2_breaks[x+1]))
    })
    
    ggplot(a) +
        scale_fill_viridis(option=virtis_option) +
        geom_tile(aes(var1, var2, fill=Prediction)) +
        xlab(var1) +
        ylab(var2) +
        theme_minimal()
}

#' run lasso
#' @export
#' @import caret

run_lasso <- function(long, variable="grouped_severity", sampling=NULL,
                   selected_features=NULL) {
    
    # convert data from long format to a matrix
    prot_matrix <- long %>% 
        dplyr::select(SampleID, GeneID, NPX, variable) %>% 
        filter(!is.na(NPX)) %>% 
        spread(GeneID, NPX) %>%
        column_to_rownames('SampleID')
    
    # select specific features for model
    if (!(is.null(selected_features))) {
        prot_matrix <- prot_matrix[,colnames(prot_matrix) %in% selected_features
                                   | colnames(prot_matrix) == variable]
    }
    
    # train model using caret
    lasso <- train(as.formula(paste0(variable, " ~ .")), 
                data = prot_matrix, 
                localImp = TRUE,
                method = "glmnet", 
                metric = "Kappa",
                proximity=TRUE,
                family="binomial",
                preProcess = c("knnImpute", "scale", "center"), 
                tuneGrid=expand.grid(alpha=1, lambda=seq(0.0001, 1, length = 20)),
                na.action = na.pass,
                trControl = trainControl(method = "repeatedcv", p = 0.7, number=10, repeats=10, sampling=sampling))

    lasso$finalModel$call$formula <- eval(lasso$call$form)
    print(lasso)
    
    return(lasso)
}
