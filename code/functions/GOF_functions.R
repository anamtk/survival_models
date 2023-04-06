

# Balanced accuracy -------------------------------------------------------


accuracy_fun <- function(predicted_df, observed_df, ID, iteration.num){
  
  df <- predicted_df %>%
    rename("Predicted_fate" = "Fate_class") %>%
    dplyr::select(-type) %>%
    filter(Iteration == iteration.num) %>%
    left_join(observed_df, by = setNames(nm = c(ID)))
  
  sumdf <- df %>%
    summarise(TP = sum((Observed_fate == 1 & Predicted_fate == 1), na.rm = T),
              FP = sum((Observed_fate == 1 & Predicted_fate == 0), na.rm = T),
              TN = sum((Observed_fate == 0 & Predicted_fate == 0), na.rm = T),
              FN = sum((Observed_fate == 0 & Predicted_fate == 1), na.rm = T),
              sensitivity = TP/(TP+FP),
              specificity = TN/(TN+FN),
              accuracy = (sensitivity+specificity)/2) %>%
    dplyr::select(sensitivity, specificity, accuracy)
  
  return(sumdf)
}



# AUC ---------------------------------------------------------------------

AUC_JAGS <- function(mod_GOF, iteration.num, resp){

  predictor_df <- as.data.frame(t(mod_GOF$sims.list$p))
  
  p1 <- predictor_df[,iteration.num]
  
  ROC <- roc(response =resp, predictor =p1)
  
  AUC <- auc(ROC)
  
  return(AUC)
}

AUC_JAGS <- function(mod_GOF, pred, iteration.num, resp){
  
  predictor_df <- as.data.frame(t(mod_GOF$sims.list$p))
  
  p1 <- predictor_df[,iteration.num]
  
  ROC <- roc(response =resp, predictor =p1)
  
  AUC <- auc(ROC)
  
  return(AUC)
}



