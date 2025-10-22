#functions for testing predictions


#CRPS (continuous rank probability score)
CPRS <- function(pred, trueobs){
  z <- as.numeric((trueobs - pred$mean) / pred$sd)
  scores <- pred$sd * (z * (2 * pnorm(z,0,1) - 1) + 
                         2 * dnorm(z,0,1) - 1/sqrt(pi))
  return(scores)
}


#IS95 (95% interval scores)
intscore <- function(x ,y, alpha = 0.05){
  hw <- -qnorm(alpha/2) * x$sd
  scores <- 2 * hw + (2/alpha) * (((x$mean - hw) - y) * (y < x$mean - hw) +
                                    (y - (x$mean +hw)) * (y > x$mean + hw))
  return(scores)
}


#coverage (coverage scores)
cvg <- function(x, y, alpha = 0.05){
  hw <- -qnorm(alpha/2) * x$sd
  scores <- y >= (x$mean - hw) & y <= (x$mean + hw)
  return(scores)
}


#TODO: finalize or move elsewhere
#press and pred R^2
predR2 <- function(model, pred, trueobs){
  
  
}


#combination of previous metrics and lm creation.
prediction_scores <- function(data_list, model_coefs, season_weeks,
                              loo = FALSE) {
  predscores <- data.frame(Week = NA, CPRS = NA, IS = NA, Coverage = NA, bic = NA )
  
  
  for(m in 1:length(season_weeks)) {
    
    temp_df <- data_list[[m]] #select correct week
    temp_df <- temp_df[-19, ] #remove 2019/2020 wildfire season
    
    
    #full data frames  
    full_pred <- temp_df[, -c(1,2)]
    full_resp <- temp_df[,2]
    
    
    #pred vars for week
    temp_pred <- full_pred[model_coefs[[m]]]
    temp_lm <- lm(full_resp ~ ., data = temp_pred)
    
    
    #for loop and setup
    preds_rf <- c()
    sd_rf <- c()
    for(i in 1:18) {
      new_data <- as.data.frame(temp_pred[i, ])
      colnames(new_data) <- colnames(temp_pred)
      
      pred_vals <- predict(temp_lm, newdata = new_data, se.fit = TRUE,
                           interval = "prediction", level = 0.95)
      preds_rf <- c(preds_rf, pred_vals$fit[1])
      sd_rf <- c(sd_rf, pred_vals$se.fit)
    }
    
    #assign the below outputs to some larger data frame
    cprs <- mean(CPRS(list(mean = preds_rf, sd = sd_rf), full_resp)) #cprs score for week 35
    is <- mean(intscore(list(mean = preds_rf, sd = sd_rf), full_resp)) #95 interval score
    cov <- mean(cvg(list(mean = preds_rf, sd = sd_rf), full_resp)) #coverage
    
    bic <- BIC(temp_lm)
    
    #add to dataframe
    predscores <- rbind(predscores, c(m, cprs, is, cov, bic))
    
  }
  
  predscores <- predscores[-1, ]
  return(predscores)
}



#add in fuzzy coverage (Paige et al 2022) Time permitting. 

