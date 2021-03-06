#' run random forests
#' @export
#' @import caret randomForestExplainer

run_rf <- function(long, variable="grouped_severity", sampling=NULL, 
                   method="rf", which_data="protein") {
    
    # convert data from long format to a matrix
    prot_matrix <- long %>% 
        dplyr::select(SampleID, GeneID, NPX, variable) %>% 
        dplyr::filter(!is.na(NPX)) %>% 
        tidyr::spread(GeneID, NPX)
    
    if (which_data %in% c("both", "clinical")) {
        
        long <- unique(dplyr::select(long, SampleID, radiology.evidence.covid, ihd, previous.vte, copd, diabetes, smoking, Sex, calc.age, macro.ethnic, cause.eskd, temp_CRP, temp_neut, temp_mono, temp_lymph, temp_WCC, temp_ddimer))
        prot_matrix <- dplyr::left_join(prot_matrix, long, by=c("SampleID"="SampleID"))
        
    }
    
    if (which_data == "clinical") {
        
        prot_matrix <- dplyr::select(prot_matrix, SampleID, radiology.evidence.covid, ihd, previous.vte, copd, grouped_severity, diabetes, smoking, Sex, calc.age, macro.ethnic, cause.eskd, temp_CRP, temp_neut, temp_mono, temp_lymph, temp_WCC, temp_ddimer)
    }
    
    prot_matrix <- tibble::column_to_rownames(prot_matrix, "SampleID")

    if (method == "rf") {
        tuneGrid <- data.frame(mtry = floor(sqrt(ncol(prot_matrix) - 1)))
        
    } else if (method == "glmnet") {
        tuneGrid <- data.frame(alpha = 0,  # tried 0.1, 0.2... 1
                               lambda = seq(0, 1, length = 1000))  # tried 0 to 10
    }
    
    # train model using caret
    rf <- train(as.formula(paste0(variable, " ~ .")), 
                data = prot_matrix, 
                localImp = TRUE,
                method = method, 
                metric = "Accuracy",
                proximity = TRUE,
                preProcess = c("knnImpute", "scale", "center"),
                na.action = na.pass,
                tuneGrid = tuneGrid,
                trControl = trainControl(method = "repeatedcv", allowParallel = TRUE, number = 4, repeats = 100, p = 0.75, sampling=sampling))
    
    rf$finalModel$call$formula <- eval(rf$call$form)
    print(rf)
    
    if (method == "rf") {
        # plot feature importance (ntrees, minimal depth)
        print(plot_min_depth_distribution(min_depth_distribution(rf$finalModel), k = 12
        ) +
            theme_minimal() +
            theme(title=element_blank(), panel.grid.major = element_blank(), text = element_text(size=11),
                  legend.position = "top") +
            guides(fill=guide_legend(title = "Number of Nodes", nrow=1,byrow=TRUE)) +
            xlab("Number of Trees"))
    }
    
    return(rf)
}

#' plot random forests
#' @export
#' @import caret randomForestExplainer
#' 
plot_rf <- function(rf, enet_imp, enet_label=0.95, mmdepth_label=2.5, 
                    return_imp_frame=FALSE) {
    
    # get feature importances
    imp_frame <- measure_importance(rf$finalModel)
    
    # adjust p values
    imp_frame$p_adj <- p.adjust(imp_frame$p_value, "fdr")
    imp_frame$p_s <- factor(ifelse(imp_frame$p_adj < 0.05, "<0.05", "NS"), levels=c("NS", "<0.05"))
    
    # use lasso importance metric
    imp_frame <- enet_imp %>%
        left_join(imp_frame, by=c("protein"="variable"))
    
    if (return_imp_frame) {
        return(imp_frame)
    }
    
    # secondary feature importance plot
    imp_plot <- ggplot(imp_frame, aes(in_x_model, mean_min_depth)) + 
        theme_pubr() +
        xlab("Elastic Net Model Inclusion Proportion") + ylab("Random Forest Mean Minimal Depth") +
        geom_point(alpha=0.75) +
        geom_text_repel(data=subset(imp_frame, in_x_model > enet_label & mean_min_depth < mmdepth_label),
                        aes(in_x_model, mean_min_depth, label = protein), size = 3, color="black")
    
    print(imp_plot)
    # print(plot_importance_ggpairs(imp_frame, c("no_of_trees", "mean_min_depth", "p_adj", "mean_min_depth", "gini_decrease")) + 
    #           theme_pubr(border = TRUE))
    
    # return(important_variables(select(imp_frame, -lasso), k = 30, measures = c("mean_min_depth", "no_of_trees")))
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
                      selected_features=NULL, lasso=TRUE) {
    
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
    
    if (lasso) {
        a <- 1
        l <- seq(0, 0.2, length = 50)
    } else {
        a <- 0
        l <- seq(0, 10, length = 50)
    }
    
    # train model using caret
    lasso <- train(as.formula(paste0(variable, " ~ .")), 
                   data = prot_matrix, 
                   localImp = TRUE,
                   method = "glmnet", 
                   metric = "Accuracy",
                   proximity=TRUE,
                   family="binomial",
                   preProcess = c("knnImpute", "scale", "center"), 
                   tuneGrid=expand.grid(alpha=a, lambda=l),
                   na.action = na.pass,
                   trControl = trainControl(method = "repeatedcv", p = 0.7, number=10, repeats=10, sampling=sampling))
    
    lasso$finalModel$call$formula <- eval(lasso$call$form)
    print(lasso)
    
    return(lasso)
}

#' get first sample for each individual
#' @export

get_first_samples <- function(long) {
    # samples c("C24", "C32", "C11") have NA for date of first symptoms, so get the first
    # samples a different way
    long$Time_From_First_Symptoms[long$SampleID %in% c("C24_211", "C11_219", "C32_188")] <- 1
    long <- long[!is.na(long$Time_From_First_Symptoms),]
    long <- long[long$Individual_ID != "C105",]  # has sample before COVID symptoms
    
    # now order the remaining samples by their time from first symptoms
    unique_SampleIDs <- unique(long$SampleID)
    unique_SampleIDs_to_remove <- vector("logical", length(unique_SampleIDs))
    
    for (ind in unique(long$Individual_ID)) {
        possible_unique_SampleIDs <- unique_SampleIDs[grepl(paste0(ind, "_"), unique_SampleIDs)]
        stopifnot(length(possible_unique_SampleIDs) > 0)
        
        if (length(possible_unique_SampleIDs) > 1) {
            SampleID_dates <- unique(data.frame(
                samp=long$SampleID[long$Individual_ID == ind],
                date=long$Time_From_First_Symptoms[long$Individual_ID == ind]))
            
            early_sample <- SampleID_dates$samp[which.min(SampleID_dates$date)]
            stopifnot(length(which(possible_unique_SampleIDs == early_sample)) == 1)
            
            possible_unique_SampleIDs <- 
                possible_unique_SampleIDs[-which(possible_unique_SampleIDs == early_sample)]
            
            stopifnot(length(possible_unique_SampleIDs) == nrow(SampleID_dates) - 1)
            unique_SampleIDs_to_remove[unique_SampleIDs %in% possible_unique_SampleIDs] <- TRUE
        }
    }
    
    independent_long <- 
        long[!(long$SampleID %in% unique_SampleIDs[unique_SampleIDs_to_remove]),]
    
    stopifnot(length(unique(independent_long$SampleID)) == length(unique(long$Individual_ID)))
    
    return(independent_long)
}
