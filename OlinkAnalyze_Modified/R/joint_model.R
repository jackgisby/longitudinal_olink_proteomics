#' Run joint models on olink data
#' @export
#' @import JMbayes

joint_model <- function(long, prot="CCL2") {
    long <- long[long$GeneID == prot,]
    
    # remove negative values of Time_From_First_Symptoms
    long <- long[long$Time_From_First_Symptoms > 0 & !is.na(long$Time_From_First_Symptoms) & long$Time_From_First_Symptoms <= 28,]
    
    # long must be sorted properly for the JM package to work
    long <- long[order(long$Individual_ID, long$Time_From_First_Symptoms),]
    
    prot <- long[long$Time_From_First_Symptoms <= max(long$Time_From_First_Symptoms[long$Fatal_Disease]),]
    prot <- prot[prot$Individual_ID != "C77",] # C77 has no measurements near date of death
    prot <- prot[prot$Individual_ID %in% names(table(prot$Individual_ID))[table(prot$Individual_ID) > 1],]  # remove single measurements
    prot$NPX <- scale(prot$NPX, center = TRUE, scale = TRUE)
    
    surv_long <- unique(data.frame(
        Individual_ID = prot$Individual_ID,
        Sex = prot$Sex,
        Age = prot$Age,
        Ethnicity = prot$Ethnicity
    ))
    
    if (any(duplicated(surv_long$Individual_ID))) {
        surv_long <- surv_long[-which(duplicated(surv_long$Individual_ID)),]
    }
    
    surv_long$is_fatal <- vector("logical", length(unique(prot$Individual_ID)))
    surv_long$time_fatal <- vector("numeric", length(unique(prot$Individual_ID)))
    
    # algorithm for finding time of event: transition to severe_critical
    # for each individual
    for (j in 1:length(unique(prot$Individual_ID))) {
        individual <- unique(prot$Individual_ID)[j]
        individual_long <- prot[prot$Individual_ID == individual,]
        individual_order <- order(individual_long$Time_From_First_Symptoms)
        
        stopifnot(length(unique(individual_long$Fatal_Disease)) == 1)
        
        if (!unique(individual_long$Fatal_Disease)) { # non-fatal
            surv_long$is_fatal[j] <- FALSE
        } else if (unique(individual_long$Fatal_Disease)) { # else: fatal
            surv_long$is_fatal[j] <- TRUE
        } else {
            stopifnot(FALSE)
        }
        
        if (surv_long$is_fatal[j]) {
            surv_long$time_fatal[j] <- as.numeric(individual_long$Time_Of_Death_From_First_Symptoms)[1]
        } else {
            surv_long$time_fatal[j] <- max(individual_long$Time_From_First_Symptoms) + 1
        }
    }
    
    surv_long <- surv_long[surv_long$Individual_ID %in% prot$Individual_ID,]
    surv_long$is_fatal <- as.numeric(surv_long$is_fatal)
    
    # fit linear mixed model
    fitLME <- lme(NPX ~ bs(Time_From_First_Symptoms, degree=2, Boundary.knots=c(0, 29)), 
                  random = list(Individual_ID = pdDiag(form = ~ bs(Time_From_First_Symptoms, degree=2, Boundary.knots=c(0, 29)))),
                  data = prot, control = lmeControl(opt='optim', msMaxIter=1000))
    
    # fit cox model
    fitSURV <- coxph(Surv(time_fatal, event=is_fatal) ~ 1, data = surv_long, x = TRUE)
    
    return(jointModelBayes(fitLME, fitSURV, timeVar = "Time_From_First_Symptoms", control=list(n.iter=20000)))
}