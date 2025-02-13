
#note: this will most likely be part of the final git repo. WRITE GOOD CODE
#TODO: produce plots for every year and region with true, pred, 95 PI (and training/test)
## -include combined plots for full data and training/testing schema
#TODO: add in residual boxplots, resid vs fitted 
#TODO: for "normal" single model validation techniques 
## -(resid v fit, qq, etc) use my partitions to make it more readable adn save space
#TODO: do prediction scores, MSE, etc to quantify preds.

#libraries


#Import Data 
#TODO: setwd to the git repo
setwd("~/CO_AUS/Aus_CO-main")

#TODO: rewrite variable_selection.rmd to only include vsurf and whatever else is used.
# importing vars as:
load( "ne_rf.rda") #(interp ne)
load( "pred_ne_rf.rda") #(pred ne)
load( "se_rf.rda") #(interp se)
load( "pred_ne_rf.rda") #(pred se)

# Data
load( "bounded_pred.rda")
load( "bounded_resp.rda")
load( "ne_data.rda")
load( "se_data.rda")

#Functions:
source("predictionTest_functions.R")

## Internal Functions (move if needed)
model_lm <- function(data_list, model_coefs, season_week){
  models_list <- list()
  for(m in 1:length(season_weeks)) {
    temp_df <- data_list[[m]] #select correct week
    
    #TODO: change this as max/length
    last_year <- length(temp_df$Seasons)
    temp_df <- temp_df[-last_year, ] #remove 2019/2020 wildfire season
    
    #full data frames  
    full_pred <- temp_df[, -c(1,2)]
    full_resp <- temp_df[,2]
    
    #pred vars for week
    temp_pred <- full_pred[model_coefs[[m]]]
    temp_lm <- lm(full_resp ~ ., data = temp_pred)
    
    models_list[[m]] <- temp_lm
  }
  return(models_list)
}


#TODO: check on sd (standard dev) output
prediction_lm <- function(model_list, data_list, season_weeks, year){
  preds_list <- c()
  upper_bound <- c()
  lower_bound <- c()
  stand_dev <- c()
  
  for (m in 1:length(season_weeks)) {
    temp_lm <- model_list[[m]]
    temp_data <- data_list[[m]]
    
    temp_coefs <- names(temp_lm$coefficients)[-1]
    model_data <- temp_data[temp_coefs]
    
    new_data <- as.data.frame(model_data[year, ])
    colnames(new_data) <- colnames(model_data)
    
    pred_vals <- predict(temp_lm, newdata = new_data, se.fit = TRUE,
                         interval = "prediction", level = 0.95)
    
    #TODO: add in sd
    preds_list <- c(preds_list, pred_vals$fit[1])
    upper_bound <- c(upper_bound, pred_vals$fit[3])
    lower_bound <- c(lower_bound, pred_vals$fit[2])
    stand_dev <- c(stand_dev, pred_vals$se.fit)
  }
  
  return(list(pred = preds_list, upper = upper_bound, 
              lower = lower_bound, std_dev = stand_dev))
}


# Other setup
season_weeks <- c(35:52, 1:14)
season_years <- c(2001:2019) #For season 2001/2002 - 2019/2020

### --- Main --- ###

# Weekly Models

#full weekly models
ne_models <- model_lm(NE_list, ne_coefs, season_weeks)
se_models <- model_lm(SE_list, se_coefs, season_weeks)


## for pred level
ne_models_pred <- model_lm(NE_list, ne_coefs_p, season_weeks)
se_models_pred <- model_lm(SE_list, ne_coefs_p, season_weeks) 

#base preds
n <- length(season_years) - 1

ne_preds <- list()
se_preds <- list()
for (i in 1:n) {
  ne_preds[[i]] <- prediction_lm(ne_models, NE_list, season_weeks, year = i)
  se_preds[[i]] <- prediction_lm(se_models, SE_list, season_weeks, year = i)
}

## for pred level

ne_pred_p <- list()
se_pred_p <- list()
for (i in 1:n) {
  ne_pred_p[[i]] <- prediction_lm(ne_models_pred, NE_list, season_weeks, year = i)
  se_pred_p[[i]] <- prediction_lm(se_models_pred, SE_list, season_weeks, year = i)
}

#set this up for weekly models
ne_preds_weekly <- list()
se_preds_weekly <- list()

for (j in 1:length(season_weeks)) {
  ne_pred_week <- c()
  ne_upper_week <- c()
  ne_lower_week <- c()
  ne_std_dev_week <- c()
  
  se_pred_week <- c()
  se_upper_week <- c()
  se_lower_week <- c()
  se_std_dev_week <- c()
  
  for (k in 1:n) {
    ne_pred_week <- c(ne_pred_week, ne_preds[[k]]$pred[j])
    ne_upper_week <- c(ne_upper_week, ne_preds[[k]]$upper[j])
    ne_lower_week <- c(ne_lower_week, ne_preds[[k]]$lower[j])
    ne_std_dev_week <- c(ne_std_dev_week, ne_preds[[k]]$std_dev[j])  
    
    se_pred_week <- c(se_pred_week, se_preds[[k]]$pred[j])
    se_upper_week <- c(se_upper_week, se_preds[[k]]$upper[j])
    se_lower_week <- c(se_lower_week, se_preds[[k]]$lower[j])
    se_std_dev_week <- c(se_std_dev_week, se_preds[[k]]$std_dev[j]) 
  }
  ne_preds_weekly[[j]] <- list(pred = ne_pred_week, upper = ne_upper_week, 
                               lower = ne_lower_week, std_dev = ne_std_dev_week)

  se_preds_weekly[[j]] <- list(pred = se_pred_week, upper = se_upper_week, 
                               lower = se_lower_week, std_dev = se_std_dev_week)    
}


#Training/Testing

remove_year <- function(df, row_year){
  df[-row_year,  ,drop = FALSE]
}

n <- length(season_years) - 1

##NE Aus train/test
ne_modeltrain <- list()
ne_predtest <- list()

for (k in 1:n) {
  NE_list_updated <- lapply(NE_list, remove_year, row_year = k)
  temp_train <- model_lm(NE_list_updated, ne_coefs, season_weeks)
  temp_test <- prediction_lm(temp_train, NE_list, season_weeks, year = k)
  
  ne_modeltrain[[k]] <- temp_train
  ne_predtest[[k]] <- temp_test
}

##SE Aus train/test
se_modeltrain <- list()
se_predtest <- list()

for (k in 1:n) {
  SE_list_updated <- lapply(SE_list, remove_year, row_year = k)
  temp_train <- model_lm(SE_list_updated, se_coefs, season_weeks)
  temp_test <- prediction_lm(temp_train, SE_list, season_weeks, year = k)
  
  se_modeltrain[[k]] <- temp_train
  se_predtest[[k]] <- temp_test
}


#yearly data (all data into a single list)
#TODO: finish SE Aus and check to make sure we didn't incorrectly re-arrange data
prediction_years <- season_years[1:18]

NE_yearly <- list()
SE_yearly <- list()

for(i in prediction_years){
  NE_year_preds <- data.frame(resp = NA, pred = NA, pred_test = NA,
                          upper = NA, lower = NA,
                          upper_test = NA, lower_test = NA)
  SE_year_preds <- data.frame(resp = NA, pred = NA, pred_test = NA,
                          upper = NA, lower = NA,
                          upper_test = NA, lower_test = NA)
  
  year <- which(prediction_years == i)
  
  #NE Aus
  for (j in 1:length(season_weeks)) {
    temp_resp <- NE_list[[j]]$NE_Aus[year]
    temp_pred <- ne_preds[[year]]$pred[j]
    temp_pred_test <- ne_predtest[[year]]$pred[j]
    
    temp_upper <- ne_preds[[year]]$upper[j]
    temp_lower <- ne_preds[[year]]$lower[j]
    
    temp_upper_test <- ne_predtest[[year]]$upper[j]
    temp_lower_test <- ne_predtest[[year]]$lower[j]
    
    NE_year_preds <- rbind(NE_year_preds, 
                           c(temp_resp, temp_pred, temp_pred_test,
                             temp_upper, temp_lower, temp_upper_test, temp_lower_test))
  } 
  
  NE_year_preds <- NE_year_preds[-1, ]
  rownames(NE_year_preds) <- season_weeks
  
  NE_yearly[[year]] <- NE_year_preds
  
  #SE Aus
  for (j in 1:length(season_weeks)) {
    temp_resp <- SE_list[[j]]$SE_Aus[year]
    temp_pred <- se_preds[[year]]$pred[j]
    temp_pred_test <- se_predtest[[year]]$pred[j]
    
    temp_upper <- se_preds[[year]]$upper[j]
    temp_lower <- se_preds[[year]]$lower[j]
    
    temp_upper_test <- se_predtest[[year]]$upper[j]
    temp_lower_test <- se_predtest[[year]]$lower[j]
    
    SE_year_preds <- rbind(SE_year_preds, 
                           c(temp_resp, temp_pred, temp_pred_test,
                             temp_upper, temp_lower, temp_upper_test, temp_lower_test))
  } 
  
  SE_year_preds <- SE_year_preds[-1, ]
  rownames(SE_year_preds) <- season_weeks
  
  SE_yearly[[year]] <- SE_year_preds
}

#save for viz elsewhere
setwd("~/CO_AUS/Aus_CO-main/Figures_Simple")
save(NE_yearly, SE_yearly, file = "simple_pred_data.rda")

#test plots:

#se aus
m <- 7 #year index
resp_test <- NE_yearly[[m]]$resp
pred_base <- NE_yearly[[m]]$pred
pred_test <- ne_pred_p[[m]]$pred
upper <- NE_yearly[[m]]$upper
lower <- NE_yearly[[m]]$lower

#prediction viz
ylim_test <- c(min(resp_test, pred_base,  pred_test), max(resp_test, pred_base,  pred_test))

png("NEAus_test2.png", width = 950, height = 500, res = 100)
par(mar = c(4, 4, 4, 12) + 0.5, xpd = TRUE)
plot(1:32, resp_test, type = "l", lwd = 2, ylim = ylim_test, xaxt = "n",
     ylab = "Anomaly (Total Column CO)", xlab = "Weeks", main = "NE Aus 2007/2008 Wildfire Season")
axis(1, at=1:32, labels = season_weeks)
lines(1:32, pred_base, lty = 2, lwd = 1.7)
lines(1:32, pred_test, lty = 4, lwd = 1.7, col = "magenta4")
legend("topright", inset = c(-0.378, 0), 
       legend = c("Actual", "Predicted (Interp Level)", "Predicted (Pred Level)"), 
       col = c("black", "black", "magenta4"), 
       lty = c(1,2,4), lwd = c(2,1.7, 1.7))
dev.off()



#residuals

#rmse
rmse_ne <- c()
rmse_se <- c()
rmse_test_ne <- c()
rmse_test_se <- c()
for (j in 1:length(season_weeks)) {
  temp_ne <- sqrt(mean((ne_models[[j]]$residuals)^2))
  temp_se <- sqrt(mean((se_models[[j]]$residuals)^2))
  
  ne_sqresid <- c()
  se_sqresid <- c()
  for (k in 1:18) {
    ne_sqresid <- c(ne_sqresid, (NE_yearly[[k]]$resp[j] - NE_yearly[[k]]$pred_test[j])^2)
    se_sqresid <- c(se_sqresid, (SE_yearly[[k]]$resp[j] - SE_yearly[[k]]$pred_test[j])^2)
  }
  rmse_ne <- c(rmse_ne, temp_ne)
  rmse_se <- c(rmse_se, temp_se)
  
  rmse_test_ne <- c(rmse_test_ne, sqrt(mean(ne_sqresid)))
  rmse_test_se <- c(rmse_test_se, sqrt(mean(se_sqresid)))
}

png("rmse.png", width = 950, height = 500, res = 100)
par(mar = c(4, 4, 4, 8) + 0.2, xpd = TRUE)
plot(1:32, rmse_ne, type = "b", pch = 16, xaxt = "n", ylim = c(0,10),
     ylab = "RMSE", xlab = "Weeks", main = "Root Mean Squared Error", 
     col = "darkorange1")
axis(1, at=1:32, labels = season_weeks)
points(1:32, rmse_se, type = "b", pch = 16, col = "magenta4")
legend("topright", inset = c(-0.128, 0), 
       legend = c("NE Aus", "SE Aus"), 
       col = c("darkorange1", "magenta4"), 
       pch = c(16, 16))
dev.off()

png("rmse_test.png", width = 950, height = 500, res = 100)
par(mar = c(4, 4, 4, 8) + 0.2, xpd = TRUE)
plot(1:32, rmse_test_ne, type = "b", pch = 16, xaxt = "n", ylim = c(0,15),
     ylab = "RMSE", xlab = "Weeks", main = "Root Mean Squared Error (Test/Train)", 
     col = "darkorange1")
axis(1, at=1:32, labels = season_weeks)
points(1:32, rmse_test_se, type = "b", pch = 16, col = "magenta4")
legend("topright", inset = c(-0.128, 0), 
       legend = c("NE Aus", "SE Aus"), 
       col = c("darkorange1", "magenta4"), 
       pch = c(16, 16))
dev.off()



##redo with new model information

#TODO: create residual histogram to see if there is systemic issues of our model
residtest <- resp_test - pred_test

#residuals of full models 
NEresid_hist <- as.data.frame(matrix(NA, nrow = 18))
SEresid_hist <- as.data.frame(matrix(NA, nrow = 18))

for (k in 1:32) {
  NEresid_hist <- cbind(NEresid_hist, ne_predlist$models[[k]]$residuals)
  SEresid_hist <- cbind(SEresid_hist, se_predlist$models[[k]]$residuals)
}
NEresid_hist <- NEresid_hist[,-1]
SEresid_hist <- SEresid_hist[,-1]

colnames(NEresid_hist) <- season_weeks[1:32]
colnames(SEresid_hist) <- season_weeks[1:32]


#resid v. fitted

multi_cols <- rainbow(8)
temp_resid <- ne_predlist$models[[9]]$residuals
temp_fit <- ne_predlist$models[[9]]$fitted.values
#get resid and fit limits across all included models (for ylim/xlim resp)

#TODO: clean up!!
par(mar = c(5, 4, 4, 8.2) + 0.5, xpd = TRUE)
plot(temp_fit, temp_resid, pch = 16, col = multi_cols[1])
abline(h = 0, lty = 2, col = "darkgrey", lwd = 2)
for (i in 10:16) {
  temp_resid <- ne_predlist$models[[i]]$residuals
  temp_fit <- ne_predlist$models[[i]]$fitted.values
  points(temp_fit, temp_resid, pch = 16, col = multi_cols[i-8])
}
#TODO: create a vector for Week #
legend("topright", inset = c(-0.32, 0),
       legend = season_weeks[9:16], col = multi_cols, pch = rep(16, length.out = 8))



#TEST resid/fit plot (looking at specific models)
test_resid <- ne_predlist$models[[3]]$residuals
test_fit <- ne_predlist$models[[3]]$fitted.values
plot(test_fit, test_resid, pch = 16)
abline(h = 0, lty = 2, col = "darkgrey", lwd = 2)


#adj R^2 and residuals
NE_adj_r2 <- c()
SE_adj_r2 <- c()
resid_test <- list()
for (k in 1:length(season_weeks)) {
  ne_model <- ne_models[[k]]
  ne_summary <- summary(ne_model)
  se_model <- se_models[[k]]
  se_summary <- summary(se_model)  
  
#  resid_test[[k]] <- temp_model$residuals
  NE_adj_r2 <- c(NE_adj_r2, ne_summary$adj.r.squared)
  SE_adj_r2 <- c(SE_adj_r2, se_summary$adj.r.squared)
}



test_summary <- summary(ne_models[[1]])
test_summary$adj.r.squared
ne_models[[1]]$residuals

## --- Model Analytics --- ##

## for resid analysis, pred scores, etc.

#residuals for weekly models
residNE1_box_df <- data.frame(values = NA, group = NA)
residSE1_box_df <- data.frame(values = NA, group = NA)
for (i in 1:16) {
  temp_vals <- ne_models[[i]]$residuals
  temp_group <- rep(season_weeks[i], length(ne_models[[i]]$residuals))
  residNE1_box_df <- rbind(residNE1_box_df, list(temp_vals, temp_group))
  
  temp_vals1 <- se_models[[i]]$residuals
  temp_group1 <- rep(season_weeks[i], length(se_models[[i]]$residuals))
  residSE1_box_df <- rbind(residSE1_box_df, list(temp_vals1, temp_group1))
}

residNE1_box_df <- residNE1_box_df[-1, ]
residNE1_box_df$group <- as.factor(residNE1_box_df$group)
residNE1_box_df$values <- as.numeric(residNE1_box_df$values)

residSE1_box_df <- residSE1_box_df[-1, ]
residSE1_box_df$group <- as.factor(residSE1_box_df$group)
residSE1_box_df$values <- as.numeric(residSE1_box_df$values)

#second half
residNE2_box_df <- data.frame(values = NA, group = NA)
residSE2_box_df <- data.frame(values = NA, group = NA)
for (i in 17:32) {
  temp_vals <- ne_models[[i]]$residuals
  temp_group <- rep(season_weeks[i], length(ne_models[[i]]$residuals))
  residNE2_box_df <- rbind(residNE2_box_df, list(temp_vals, temp_group))
  
  temp_vals1 <- se_models[[i]]$residuals
  temp_group1 <- rep(season_weeks[i], length(se_models[[i]]$residuals))
  residSE2_box_df <- rbind(residSE2_box_df, list(temp_vals1, temp_group1))
}

residNE2_box_df <- residNE2_box_df[-1, ]
residNE2_box_df$group <- as.factor(residNE2_box_df$group)
residNE2_box_df$values <- as.numeric(residNE2_box_df$values)

residSE2_box_df <- residSE2_box_df[-1, ]
residSE2_box_df$group <- as.factor(residSE2_box_df$group)
residSE2_box_df$values <- as.numeric(residSE2_box_df$values)

resid_min <- min(residNE1_box_df$values, residNE2_box_df$values, 
                 residSE1_box_df$values, residSE2_box_df$values)

resid_max <- max(residNE1_box_df$values, residNE2_box_df$values, 
                 residSE1_box_df$values, residSE2_box_df$values)

#resid_lim <- c(resid_min, resid_max)

resisNE_base <- rbind(residNE1_box_df ,residNE2_box_df)
resisSE_base <- rbind(residSE1_box_df ,residSE2_box_df)

resid_lim <- c(-20, 20)

#resid boxplots
setwd("~/CO_AUS/Aus_CO-main/Figures_Simple")


png("NEresid_base.png", width = 800, height = 600, res = 100)
par(mar = c(8, 4, 4, 2) + 0.1)
boxplot(values ~ group, data = resisNE_base, ylim = resid_lim, ylab = "Residuals", xlab = "",
        main = "NE Aus Model Residuals",axes = FALSE, pch = 20)
box()
axis(2)
axis(1, at = 1:32, labels = c(paste0("Week ", season_weeks)),
     las = 3)
abline(h = 0, lty = 2)
dev.off()

png("SEresid_base.png", width = 800, height = 600, res = 100)
par(mar = c(8, 4, 4, 2) + 0.1)
boxplot(values ~ group, data = resisSE_base, ylim = resid_lim, ylab = "Residuals", xlab = "",
        main = "SE Aus Model Residuals",axes = FALSE, pch = 20)
box()
axis(2)
axis(1, at = 1:32, labels = c(paste0("Week ", season_weeks)),
     las = 3)
abline(h = 0, lty = 2)
dev.off()



png("NEresid_box1.png", width = 800, height = 600, res = 100)
par(mar = c(8, 4, 4, 2) + 0.1)
boxplot(values ~ group, data = residNE1_box_df, ylim = resid_lim, ylab = "Residuals", xlab = "",
        main = "NE Aus Model Residuals (Weeks 35 to 50)",axes = FALSE, pch = 20)
box()
axis(2)
axis(1, at = 1:16, labels = c(paste0("Week ", season_weeks[1:16])),
     las = 3)
abline(h = 0, lty = 2)
dev.off()


png("SEresid_box1.png", width = 800, height = 600, res = 100)
par(mar = c(8, 4, 4, 2) + 0.1)
boxplot(values ~ group, data = residSE1_box_df, ylim = resid_lim, ylab = "Residuals", xlab = "",
        main = "SE Aus Model Residuals (Weeks 35 to 50)",axes = FALSE, pch = 20)
box()
axis(2)
axis(1, at = 1:16, labels = c(paste0("Week ", season_weeks[1:16])),
     las = 3)
abline(h = 0, lty = 2)
dev.off()


png("NEresid_box2.png", width = 800, height = 600, res = 100)
par(mar = c(8, 4, 4, 2) + 0.1)
boxplot(values ~ group, data = residNE2_box_df, ylim = resid_lim, ylab = "Residuals", xlab = "",
        main = "NE Aus Model Residuals (Weeks 51, 52, 1 to 14)",axes = FALSE, pch = 20)
box()
axis(2)
axis(1, at = 1:16, labels = c(paste0("Week ", season_weeks[17:32])),
     las = 3)
abline(h = 0, lty = 2)
dev.off()


png("SEresid_box2.png", width = 800, height = 600, res = 100)
par(mar = c(8, 4, 4, 2) + 0.1)
boxplot(values ~ group, data = residSE2_box_df, ylim = resid_lim, ylab = "Residuals", xlab = "",
        main = "SE Aus Model Residuals (Weeks 51, 52, 1 to 14)",axes = FALSE, pch = 20)
box()
axis(2)
axis(1, at = 1:16, labels = c(paste0("Week ", season_weeks[17:32])),
     las = 3)
abline(h = 0, lty = 2)
dev.off()

# Prediction Scores
#Using CPRS(), intscore(), and cvg()
#we need predictions and standard dev


#continuous probability ranked score (CPRS)
cprs_df <- data.frame(NE_Aus = NA, SE_Aus = NA)
for (j in 1:length(season_weeks)) {
  ne_resp <- NE_list[[j]][-19, 2]
  se_resp <- SE_list[[j]][-19, 2]
  
  ne_cprs <- mean(CPRS(list(mean = ne_preds_weekly[[j]]$pred, sd = ne_preds_weekly[[j]]$std_dev), ne_resp)) 
  se_cprs <- mean(CPRS(list(mean = se_preds_weekly[[j]]$pred, sd = se_preds_weekly[[j]]$std_dev), se_resp))
  
  cprs_df <- rbind(cprs_df, c(ne_cprs, se_cprs))
}
cprs_df <- cprs_df[-1, ]


png("CRPS.png", width = 950, height = 500, res = 100)
par(mar = c(4, 4, 4, 8) + 0.2, xpd = TRUE)
plot(1:32, cprs_df$NE_Aus, type = "b", pch = 16, xaxt = "n", ylim = c(0,6),
     ylab = "CRPS", xlab = "Weeks", main = "Continuous Ranked Probability Score", 
     col = "darkorange1")
axis(1, at=1:32, labels = season_weeks)
points(1:32, cprs_df$SE_Aus, type = "b", pch = 16, col = "magenta4")
legend("topright", inset = c(-0.128, 0), 
       legend = c("NE Aus", "SE Aus"), 
       col = c("darkorange1", "magenta4"), 
       pch = c(16, 16))
dev.off()

#interval score (is)
is_df <- data.frame(NE_Aus = NA, SE_Aus = NA)
for (j in 1:length(season_weeks)) {
  ne_resp <- NE_list[[j]][-19, 2]
  se_resp <- SE_list[[j]][-19, 2]
  
  ne_is <- mean(intscore(list(mean = ne_preds_weekly[[j]]$pred, sd = ne_preds_weekly[[j]]$std_dev), ne_resp)) 
  se_is <- mean(intscore(list(mean = se_preds_weekly[[j]]$pred, sd = se_preds_weekly[[j]]$std_dev), se_resp))
  
  is_df <- rbind(is_df, c(ne_is, se_is))
}
is_df <- is_df[-1, ]


png("intscore.png", width = 950, height = 500, res = 100)
par(mar = c(4, 4, 4, 8) + 0.2, xpd = TRUE)
plot(1:32, is_df$NE_Aus, type = "b", pch = 16, xaxt = "n", ylim = c(0,100),
     ylab = "Interval Score", xlab = "Weeks", main = "95% Interval Score",
     col = "darkorange1")
axis(1, at=1:32, labels = season_weeks)
points(1:32, is_df$SE_Aus, type = "b", pch = 16, col = "magenta4")
legend("topright", inset = c(-0.128, 0), 
       legend = c("NE Aus", "SE Aus"), 
       col = c("darkorange1", "magenta4"), 
       pch = c(16, 16))
dev.off()


#coverage (cvg)
cvg_df <- data.frame(NE_Aus = NA, SE_Aus = NA)
for (j in 1:length(season_weeks)) {
  ne_resp <- NE_list[[j]][-19, 2]
  se_resp <- SE_list[[j]][-19, 2]
  
  ne_cvg <- mean(cvg(list(mean = ne_preds_weekly[[j]]$pred, sd = ne_preds_weekly[[j]]$std_dev), ne_resp)) 
  se_cvg <- mean(cvg(list(mean = se_preds_weekly[[j]]$pred, sd = se_preds_weekly[[j]]$std_dev), se_resp))
  
  cvg_df <- rbind(cvg_df, c(ne_cvg, se_cvg))
}
cvg_df <- cvg_df[-1, ]


#coverage plot
png("coverage.png", width = 950, height = 500, res = 100)
par(mar = c(4, 4, 4, 8) + 0.2, xpd = TRUE)
plot(1:32, cvg_df$NE_Aus, type = "b", pch = 16, xaxt = "n", ylim = c(0,1),
     ylab = "Coverage", xlab = "Weeks", main = "95% Coverage",
     col = "darkorange1")
axis(1, at=1:32, labels = season_weeks)
points(1:32, cvg_df$SE_Aus, type = "b", pch = 16, col = "magenta4")
clip(0, 33, 0, 1) 
abline(h = 0.95, lty = 2, lwd = 1.5, col = "darkgray")
par(xpd = NA)
legend("topright", inset = c(-0.128, 0), 
       legend = c("NE Aus", "SE Aus"), 
       col = c("darkorange1", "magenta4"), 
       pch = c(16, 16))
dev.off()




