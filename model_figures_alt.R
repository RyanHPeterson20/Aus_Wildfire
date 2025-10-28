#figures from models and model validation

#notes: 
## create figures associated with models

#TODO: create sections for figures (e.g. numbered)
#1. 
# . Predictors and Model Fits
# . Predictions and Intervals
# . Adj. R^2 
# . K-fold CV plots
# . RMSE plots (see line ..)
# . CRPS plots
#n. Interval Score Plots


#libraries
suppressMessages( library(grid)) #gridlines between plots
suppressMessages( library(scales)) #for adjusting opacity

suppressMessages( library(RColorBrewer)) #colorRampPalette
suppressMessages( library( fields)) #for set.panel() and others 

#data/model imports
setwd("~/CO_AUS/Aus_CO-main/Interactions_New")
load("validation_refits.rda") #refits and validation (BIC)
load("validation_refitsEBIC.rda") #refits and validation (eBIC)
load("validation_kfold.rda") #kfold cv for both BIC and eBIC
load("validation_predR2.rda") #prediction R2 (R-squared/P-squared)
load("validation_refits_alt.rda") #refits and validation from alternative fits 
load("validation_refits_noOLR.rda") #temporary refits for non-OLR preds

setwd("~/CO_AUS/Aus_CO-main/Interactions")
source("group_functionsNew.R") #new groupings


## --- setup --- ##
#normal/base 
#season years/weeks
season_weeks <- c(35:52, 1:14)
season_years <- unique(bounded_resp_df$year)

seasons <- c()
for (i in 1:(length(season_years)-1)) {
  temp_season <- paste0(season_years[i], "-", season_years[i+1])
  #print(temp_season)  
  seasons <- c(seasons, temp_season)
}
rm(i, temp_season)

#center response data
NEbase_matrix <- scale(resp_matrix[ ,1:32], center = TRUE, scale = FALSE)
SEbase_matrix <- scale(resp_matrix[ ,33:64], center = TRUE, scale = FALSE)

#TODO: update below into a helper function
#Below includes removal of weeks 35-37 (adjust grouping functions if we need to change this.)
#NEAus
NEAus_1 <- NEbase_matrix[ ,4:12] 
NEAus_2 <- NEbase_matrix[ ,13:17]
NEAus_3 <- NEbase_matrix[ ,18:32]

#NEAus
SEAus_1 <- SEbase_matrix[ ,4:16] 
SEAus_2 <- SEbase_matrix[ ,17:20] 
SEAus_3 <- SEbase_matrix[ ,21:32]

NEAus_mat <- list(NEAus_1, NEAus_2, NEAus_3)
SEAus_mat <- list(SEAus_1, SEAus_2, SEAus_3)

rm(NEAus_1, NEAus_2, NEAus_3, SEAus_1, SEAus_2, SEAus_3)

#full model (uses `group_functionsNew.R`)
NEpreds_new <- NElag_3group(NE_laglist = NE_iodlag, j = 1:19)
NEresp_new <- NEresp_3group(NEAus_mat = NEAus_mat, j = 1:19)

SEpreds_new <- SElag_3group(SE_laglist = SE_iodlag, j = 1:19)
SEresp_new <- SEresp_3group(SEAus_mat = SEAus_mat, j = 1:19)

#TODO: move model fits up here, if needed (check later)

#TODO: get the changes in the number of terms for eBIC models
new.season.weeks <- season_weeks[-c(1:3)]


# . (n.) Predictions and Intervals
#get predictions and associated 95\% intervals

#TODO: finalize plots with legends and output as .png

#NE Aus predictions
NEpreds.alt1 <- NEvalid.alt1[[4]]




#SE Aus predictions
SEpreds.alt1 <- SEvalid.alt1[[4]]
SEpreds.alt2 <- SEvalid.alt2[[4]] #noOLR

#2006-2007 Season
temp.2006.preds <- SEpreds.alt1$`2006-2007`
plot(1:29, temp.2006.preds$true, type = "l", ylim = range(temp.2006.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2006.preds$base.fit, lty = 2)
#lines(1:29,  temp.2006.preds$const.lwr, lty = 2, col = "royalblue3")
#lines(1:29,  temp.2006.preds$const.upr, lty = 2, col = "firebrick3")
lines(1:29, temp.2006.preds$const.fit, lty = 4, lwd = 1.82)
#lines(1:29,  temp.2019.preds$const.lwr, lty = 4, lwd = 1, col = "royalblue3")
#lines(1:29,  temp.2019.preds$const.upr, lty = 4, lwd = 1, col = "firebrick3")
lines(1:29, temp.2006.preds$vary.fit, lty = 6, lwd = 1.42, col = "darkmagenta")
abline(v = c(13.5, 17.5), lty = 3, lwd = 0.75)
title("SE Aus : 2006-2007 Season", adj = 0)
abline(h=0, lty = 3)


temp.2006.preds <- SEpreds.alt2$`2006-2007`
plot(1:29, temp.2006.preds$true, type = "l", ylim = range(temp.2006.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2006.preds$base.fit, lty = 2)
#lines(1:29,  temp.2006.preds$const.lwr, lty = 2, col = "royalblue3")
#lines(1:29,  temp.2006.preds$const.upr, lty = 2, col = "firebrick3")
lines(1:29, temp.2006.preds$const.fit, lty = 4, lwd = 1.82)
#lines(1:29,  temp.2019.preds$const.lwr, lty = 4, lwd = 1, col = "royalblue3")
#lines(1:29,  temp.2019.preds$const.upr, lty = 4, lwd = 1, col = "firebrick3")
lines(1:29, temp.2006.preds$vary.fit, lty = 6, lwd = 1.22, col = "darkmagenta")
abline(v = c(13.5, 17.5), lty = 3, lwd = 0.75)
title("SE Aus : 2006-2007 Season", adj = 0)
abline(h=0, lty = 3)


par(mar = c(5, 5, 5, 4.3))
temp.2019.preds <- SEpreds.alt1$`2019-2020`
plot(1:29, temp.2019.preds$true, type = "l", ylim = range(temp.2019.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE, lwd = 1.52)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.95)
axis(2)  
lines(1:29, temp.2019.preds$base.fit, lty = 2, lwd = 1.82)
#lines(1:29,  temp.2019.preds$base.lwr, lty = 2, lwd = 1, col = "royalblue3")
#lines(1:29,  temp.2019.preds$base.upr, lty = 2, lwd = 1, col = "firebrick3")
lines(1:29, temp.2019.preds$const.fit, lty = 4, lwd = 1.82)
#lines(1:29,  temp.2019.preds$const.lwr, lty = 4, lwd = 1, col = "royalblue3")
#lines(1:29,  temp.2019.preds$const.upr, lty = 4, lwd = 1, col = "firebrick3")
lines(1:29, temp.2019.preds$vary.fit, lty = 6, lwd = 1.82, col = "darkmagenta")
abline(h=0, lty = 3)

par(mar = c(5, 5, 5, 4.3))
temp.2019.preds <- SEpreds.alt2$`2019-2020`
plot(1:29, temp.2019.preds$true, type = "l", ylim = range(temp.2019.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE, lwd = 1.52)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.95)
axis(2)  
lines(1:29, temp.2019.preds$base.fit, lty = 2, lwd = 1.82)
#lines(1:29,  temp.2019.preds$base.lwr, lty = 2, lwd = 1, col = "royalblue3")
#lines(1:29,  temp.2019.preds$base.upr, lty = 2, lwd = 1, col = "firebrick3")
lines(1:29, temp.2019.preds$const.fit, lty = 4, lwd = 1.82)
#lines(1:29,  temp.2019.preds$const.lwr, lty = 4, lwd = 1, col = "royalblue3")
#lines(1:29,  temp.2019.preds$const.upr, lty = 4, lwd = 1, col = "firebrick3")
lines(1:29, temp.2019.preds$vary.fit, lty = 6, lwd = 1.82, col = "darkmagenta")
abline(h=0, lty = 3)


