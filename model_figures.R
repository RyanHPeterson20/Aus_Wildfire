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

## --- Main --- ##


# . (n.) Predictors and Model Fits
## Not sure what I hope to accomplish here, these are very complex models
#Plot relationship between single predictor and response, holding all other vars const.
#Start with only significant predictors
#TODO: add interaction terms as x,y,z plot (fit surface according to interaction and each individual term)


#NE Aus Group 1
summary(NEmodels[[1]])
NE.y1 <- NEresp_new[[1]]
NE.X1 <- NEpreds_new[[1]]

#NE Aus Group 2

#NE Aus Group 3


#SE Aus Group 1
summary(SEmodels[[1]])
SE.y1 <- SEresp_new[[1]]
SE.X1 <- SEpreds_new[[1]]

#wtio_lag13 (sig. w/ 0.016033)
plot(SE.X1$wtio_lag13, SE.y1, pch = 16)
#etio_lag31 (sig. w/ 0.035191 )
plot(SE.X1$etio_lag31, SE.y1, pch = 16)
#aao_lag24 (sig. w/ 1.3e-06)
plot(SE.X1$aao_lag24, SE.y1, pch = 16)
#aao_lag29 (sig. w/ )
plot(SE.X1$aao_lag29, SE.y1, pch = 16)
#aao_lag35 (sig. w/ )
plot(SE.X1$aao_lag45, SE.y1, pch = 16)
#SEolr_lag28 (sig. w/ )
plot(SE.X1$SEolr_lag28, SE.y1, pch = 16)

#SE Aus Group 2
summary(SEmodels[[2]])
SE.y2 <- SEresp_new[[2]]
SE.X2 <- SEpreds_new[[2]]

#nino_lag40 (sig. w/ )
plot(SE.X2$nino_lag40, SE.y2, pch = 16)
#etio_lag7 (sig. w/ )
plot(SE.X2$etio_lag7, SE.y2, pch = 16)

#SE Aus Group 3
summary(SEmodels[[3]])
SE.y3 <- SEresp_new[[3]]
SE.X3 <- SEpreds_new[[3]]

#wtio_lag1 (sig. w/ )
plot(SE.X3$wtio_lag1, SE.y3, pch = 16)
#etio_lag15 (sig. w/ )
plot(SE.X3$etio_lag15, SE.y3, pch = 16)
#tsa_lag22 (sig. w/ )
plot(SE.X3$tsa_lag22, SE.y3, pch = 16)
#SEolr_lag1 (sig. w/ )
plot(SE.X3$SEolr_lag1, SE.y3, pch = 16)

# . (n.) Predictions and Intervals
#get predictions and associated 95\% intervals

#TODO: finalize plots with legends and output as .png

#NE Aus predictions
NEpreds <- NEvalid[[4]]

#2001-2002 Season
temp.2001.preds <- NEpreds$`2001-2002`
plot(1:29, temp.2001.preds$true, type = "l", ylim = range(temp.2001.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2001.preds$const.fit, lty = 2)
lines(1:29,  temp.2001.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2001.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(9.5, 14.5), lty = 3, lwd = 0.75)
title("NE Aus : 2001-2002 Season", adj = 0)

#2002-2003 Season
temp.2002.preds <- NEpreds$`2002-2003`
plot(1:29, temp.2002.preds$true, type = "l", ylim = range(temp.2002.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2002.preds$const.fit, lty = 2)
lines(1:29,  temp.2002.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2002.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(9.5, 14.5), lty = 3, lwd = 0.75)
title("NE Aus : 2002-2003 Season", adj = 0)

#2003-2004 Season
temp.2003.preds <- NEpreds$`2003-2004`
plot(1:29, temp.2003.preds$true, type = "l", ylim = range(temp.2003.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2003.preds$const.fit, lty = 2)
lines(1:29,  temp.2003.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2003.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(9.5, 14.5), lty = 3, lwd = 0.75)
title("NE Aus : 2003-2004 Season", adj = 0)

#2004-2005 Season
temp.2004.preds <- NEpreds$`2004-2005`
plot(1:29, temp.2004.preds$true, type = "l", ylim = range(temp.2004.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2004.preds$const.fit, lty = 2)
lines(1:29,  temp.2004.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2004.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(9.5, 14.5), lty = 3, lwd = 0.75)
title("NE Aus : 2004-2005 Season", adj = 0)

#2005-2006 Season
temp.2005.preds <- NEpreds$`2005-2006`
plot(1:29, temp.2005.preds$true, type = "l", ylim = range(temp.2005.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2005.preds$const.fit, lty = 2)
lines(1:29,  temp.2005.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2005.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(9.5, 14.5), lty = 3, lwd = 0.75)
title("NE Aus : 2005-2006 Season", adj = 0)

#2006-2007 Season
temp.2006.preds <- NEpreds$`2006-2007`
plot(1:29, temp.2006.preds$true, type = "l", ylim = range(temp.2006.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2006.preds$const.fit, lty = 2)
lines(1:29,  temp.2006.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2006.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(9.5, 14.5), lty = 3, lwd = 0.75)
title("NE Aus : 2006-2007 Season", adj = 0)

#2007-2008 Season
temp.2007.preds <- NEpreds$`2007-2008`
plot(1:29, temp.2007.preds$true, type = "l", ylim = range(temp.2007.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2007.preds$const.fit, lty = 2)
lines(1:29,  temp.2007.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2007.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(9.5, 14.5), lty = 3, lwd = 0.75)
title("NE Aus : 2007-2008 Season", adj = 0)

#2008-2009 Season
temp.2008.preds <- NEpreds$`2008-2009`
plot(1:29, temp.2008.preds$true, type = "l", ylim = range(temp.2008.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2008.preds$const.fit, lty = 2)
lines(1:29,  temp.2008.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2008.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(9.5, 14.5), lty = 3, lwd = 0.75)
title("NE Aus : 2008-2009 Season", adj = 0)

#2009-2010 Season
temp.2009.preds <- NEpreds$`2009-2010`
plot(1:29, temp.2009.preds$true, type = "l", ylim = range(temp.2009.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2009.preds$const.fit, lty = 2)
lines(1:29,  temp.2009.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2009.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(9.5, 14.5), lty = 3, lwd = 0.75)
title("NE Aus : 2009-2010 Season", adj = 0)

#2010-2011 Season
temp.2010.preds <- NEpreds$`2010-2011`
plot(1:29, temp.2010.preds$true, type = "l", ylim = range(temp.2010.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2010.preds$const.fit, lty = 2)
lines(1:29,  temp.2010.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2010.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(9.5, 14.5), lty = 3, lwd = 0.75)
title("NE Aus : 2010-2011 Season", adj = 0)

#2011-2012 Season
temp.2011.preds <- NEpreds$`2011-2012`
plot(1:29, temp.2011.preds$true, type = "l", ylim = range(temp.2011.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2011.preds$const.fit, lty = 2)
lines(1:29,  temp.2011.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2011.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(9.5, 14.5), lty = 3, lwd = 0.75)
title("NE Aus : 2011-2012 Season", adj = 0)

#2012-2013 Season
temp.2012.preds <- NEpreds$`2012-2013`
plot(1:29, temp.2012.preds$true, type = "l", ylim = range(temp.2012.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2012.preds$const.fit, lty = 2)
lines(1:29,  temp.2012.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2012.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(9.5, 14.5), lty = 3, lwd = 0.75)
title("NE Aus : 2012-2013 Season", adj = 0)

#2013-2014 Season
temp.2013.preds <- NEpreds$`2013-2014`
plot(1:29, temp.2013.preds$true, type = "l", ylim = range(temp.2013.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2013.preds$const.fit, lty = 2)
lines(1:29,  temp.2013.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2013.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(9.5, 14.5), lty = 3, lwd = 0.75)
title("NE Aus : 2013-2014 Season", adj = 0)

#2014-2015 Season
temp.2014.preds <- NEpreds$`2014-2015`
plot(1:29, temp.2014.preds$true, type = "l", ylim = range(temp.2014.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2014.preds$const.fit, lty = 2)
lines(1:29,  temp.2014.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2014.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(9.5, 14.5), lty = 3, lwd = 0.75)
title("NE Aus : 2014-2015 Season", adj = 0)

#2015-2016 Season
temp.2015.preds <- NEpreds$`2015-2016`
plot(1:29, temp.2015.preds$true, type = "l", ylim = range(temp.2015.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2015.preds$const.fit, lty = 2)
lines(1:29,  temp.2015.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2015.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(9.5, 14.5), lty = 3, lwd = 0.75)
title("NE Aus : 2015-2016 Season", adj = 0)

#2016-2017 Season
temp.2016.preds <- NEpreds$`2016-2017`
plot(1:29, temp.2016.preds$true, type = "l", ylim = range(temp.2016.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2016.preds$const.fit, lty = 2)
lines(1:29,  temp.2016.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2016.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(9.5, 14.5), lty = 3, lwd = 0.75)
title("NE Aus : 2016-2017 Season", adj = 0)

#2017-2018 Season
temp.2017.preds <- NEpreds$`2017-2018`
plot(1:29, temp.2017.preds$true, type = "l", ylim = range(temp.2017.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2017.preds$const.fit, lty = 2)
lines(1:29,  temp.2017.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2017.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(9.5, 14.5), lty = 3, lwd = 0.75)
title("NE Aus : 2017-2018 Season", adj = 0)

#2018-2019 Season
temp.2018.preds <- NEpreds$`2018-2019`
plot(1:29, temp.2018.preds$true, type = "l", ylim = range(temp.2018.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2018.preds$const.fit, lty = 2)
lines(1:29,  temp.2018.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2018.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(9.5, 14.5), lty = 3, lwd = 0.75)
title("NE Aus : 2018-2019 Season", adj = 0)

#2019-2020 Season
temp.2019.preds <- NEpreds$`2019-2020`
plot(1:29, temp.2019.preds$true, type = "l", ylim = range(temp.2019.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2019.preds$const.fit, lty = 2)
lines(1:29,  temp.2019.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2019.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(9.5, 14.5), lty = 3, lwd = 0.75)
title("NE Aus : 2019-2020 Season", adj = 0)

#SE Aus predictions
SEpreds <- SEvalid[[4]]

#2001-2002 Season
temp.2001.preds <- SEpreds$`2001-2002`
plot(1:29, temp.2001.preds$true, type = "l", ylim = range(temp.2001.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2001.preds$const.fit, lty = 2)
lines(1:29,  temp.2001.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2001.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(13.5, 17.5), lty = 3, lwd = 0.75)
title("SE Aus : 2001-2002 Season", adj = 0)

#2002-2003 Season
temp.2002.preds <- SEpreds$`2002-2003`
plot(1:29, temp.2002.preds$true, type = "l", ylim = range(temp.2002.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2002.preds$const.fit, lty = 2)
lines(1:29,  temp.2002.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2002.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(13.5, 17.5), lty = 3, lwd = 0.75)
title("SE Aus : 2002-2003 Season", adj = 0)

#2003-2004 Season
temp.2003.preds <- SEpreds$`2003-2004`
plot(1:29, temp.2003.preds$true, type = "l", ylim = range(temp.2003.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2003.preds$const.fit, lty = 2)
lines(1:29,  temp.2003.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2003.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(13.5, 17.5), lty = 3, lwd = 0.75)
title("SE Aus : 2003-2004 Season", adj = 0)

#2004-2005 Season
temp.2004.preds <- SEpreds$`2004-2005`
plot(1:29, temp.2004.preds$true, type = "l", ylim = range(temp.2004.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2004.preds$const.fit, lty = 2)
lines(1:29,  temp.2004.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2004.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(13.5, 17.5), lty = 3, lwd = 0.75)
title("SE Aus : 2004-2005 Season", adj = 0)

#2005-2006 Season
temp.2005.preds <- SEpreds$`2005-2006`
plot(1:29, temp.2005.preds$true, type = "l", ylim = range(temp.2005.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2005.preds$const.fit, lty = 2)
lines(1:29,  temp.2005.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2005.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(13.5, 17.5), lty = 3, lwd = 0.75)
title("SE Aus : 2005-2006 Season", adj = 0)

#2006-2007 Season
temp.2006.preds <- SEpreds$`2006-2007`
plot(1:29, temp.2006.preds$true, type = "l", ylim = range(temp.2006.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2006.preds$const.fit, lty = 2)
lines(1:29,  temp.2006.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2006.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(13.5, 17.5), lty = 3, lwd = 0.75)
title("SE Aus : 2006-2007 Season", adj = 0)

#2007-2008 Season
temp.2007.preds <- SEpreds$`2007-2008`
plot(1:29, temp.2007.preds$true, type = "l", ylim = range(temp.2007.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2007.preds$const.fit, lty = 2)
lines(1:29,  temp.2007.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2007.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(13.5, 17.5), lty = 3, lwd = 0.75)
title("SE Aus : 2007-2008 Season", adj = 0)

#2008-2009 Season
temp.2008.preds <- SEpreds$`2008-2009`
plot(1:29, temp.2008.preds$true, type = "l", ylim = range(temp.2008.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2008.preds$const.fit, lty = 2)
lines(1:29,  temp.2008.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2008.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(13.5, 17.5), lty = 3, lwd = 0.75)
title("SE Aus : 2008-2009 Season", adj = 0)

#2009-2010 Season
temp.2009.preds <- SEpreds$`2009-2010`
plot(1:29, temp.2009.preds$true, type = "l", ylim = range(temp.2009.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2009.preds$const.fit, lty = 2)
lines(1:29,  temp.2009.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2009.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(13.5, 17.5), lty = 3, lwd = 0.75)
title("SE Aus : 2009-2010 Season", adj = 0)

#2010-2011 Season
temp.2010.preds <- SEpreds$`2010-2011`
plot(1:29, temp.2010.preds$true, type = "l", ylim = range(temp.2010.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2010.preds$const.fit, lty = 2)
lines(1:29,  temp.2010.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2010.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(13.5, 17.5), lty = 3, lwd = 0.75)
title("SE Aus : 2010-2011 Season", adj = 0)

#2011-2012 Season
temp.2011.preds <- SEpreds$`2011-2012`
plot(1:29, temp.2011.preds$true, type = "l", ylim = range(temp.2011.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2011.preds$const.fit, lty = 2)
lines(1:29,  temp.2011.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2011.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(13.5, 17.5), lty = 3, lwd = 0.75)
title("SE Aus : 2011-2012 Season", adj = 0)

#2012-2013 Season
temp.2012.preds <- SEpreds$`2012-2013`
plot(1:29, temp.2012.preds$true, type = "l", ylim = range(temp.2012.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2012.preds$const.fit, lty = 2)
lines(1:29,  temp.2012.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2012.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(13.5, 17.5), lty = 3, lwd = 0.75)
title("SE Aus : 2012-2013 Season", adj = 0)

#2013-2014 Season
temp.2013.preds <- SEpreds$`2013-2014`
plot(1:29, temp.2013.preds$true, type = "l", ylim = range(temp.2013.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2013.preds$const.fit, lty = 2)
lines(1:29,  temp.2013.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2013.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(13.5, 17.5), lty = 3, lwd = 0.75)
title("SE Aus : 2013-2014 Season", adj = 0)

#2014-2015 Season
temp.2014.preds <- SEpreds$`2014-2015`
plot(1:29, temp.2014.preds$true, type = "l", ylim = range(temp.2014.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2014.preds$const.fit, lty = 2)
lines(1:29,  temp.2014.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2014.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(13.5, 17.5), lty = 3, lwd = 0.75)
title("SE Aus : 2014-2015 Season", adj = 0)

#2015-2016 Season
temp.2015.preds <- SEpreds$`2015-2016`
plot(1:29, temp.2015.preds$true, type = "l", ylim = range(temp.2015.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2015.preds$const.fit, lty = 2)
lines(1:29,  temp.2015.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2015.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(13.5, 17.5), lty = 3, lwd = 0.75)
title("SE Aus : 2015-2016 Season", adj = 0)

#2016-2017 Season
temp.2016.preds <- SEpreds$`2016-2017`
plot(1:29, temp.2016.preds$true, type = "l", ylim = range(temp.2016.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2016.preds$const.fit, lty = 2)
lines(1:29,  temp.2016.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2016.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(13.5, 17.5), lty = 3, lwd = 0.75)
title("SE Aus : 2016-2017 Season", adj = 0)

#2017-2018 Season
temp.2017.preds <- SEpreds$`2017-2018`
plot(1:29, temp.2017.preds$true, type = "l", ylim = range(temp.2017.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2017.preds$const.fit, lty = 2)
lines(1:29,  temp.2017.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2017.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(13.5, 17.5), lty = 3, lwd = 0.75)
title("SE Aus : 2017-2018 Season", adj = 0)

#2018-2019 Season
temp.2018.preds <- SEpreds$`2018-2019`
plot(1:29, temp.2018.preds$true, type = "l", ylim = range(temp.2018.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2018.preds$const.fit, lty = 2)
lines(1:29,  temp.2018.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2018.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(13.5, 17.5), lty = 3, lwd = 0.75)
title("SE Aus : 2018-2019 Season", adj = 0)

#2019-2020 Season
temp.2019.preds <- SEpreds$`2019-2020`
plot(1:29, temp.2019.preds$true, type = "l", ylim = range(temp.2019.preds),
     xlab = "Week", ylab = "CO Anomaly", axes = FALSE)
box()
axis(1, labels = new.season.weeks, at = 1:29, cex.axis = 0.75)
axis(2)  
lines(1:29, temp.2019.preds$const.fit, lty = 2)
lines(1:29,  temp.2019.preds$const.lwr, lty = 2, col = "royalblue3")
lines(1:29,  temp.2019.preds$const.upr, lty = 2, col = "firebrick3")
abline(v = c(13.5, 17.5), lty = 3, lwd = 0.75)
title("SE Aus : 2019-2020 Season", adj = 0)


#TODO: repeat with the points plots for the entire study period.
## See either of the two preceding papers for examples.


# . (n.) adj R^2 (in-sample)
#TODO: finalize the below figures

#Get all seasons
#NE Aus 
NE.cons.adjR2 <- NULL
NE.vary.adjR2 <- NULL
for (j in 1:3) {
  cons.adjR2 <- numeric(19)
  vary.adjR2 <- numeric(19)
  for (i in 1:19) {
    cons.lm.temp <- NErefit.new[[2]][[i]]
    vary.lm.temp <- NErefit.new[[3]][[i]]

    cons.adjR2[i] <- summary(cons.lm.temp[[j]])$adj.r.squared
    vary.adjR2[i] <- summary(vary.lm.temp[[j]])$adj.r.squared
  }
  NE.cons.adjR2[[j]] <- cons.adjR2
  NE.vary.adjR2[[j]] <- vary.adjR2
}

#NE Aus group 1
plot(1:19, NE.cons.adjR2[[1]], pch = 16, col = "firebrick", ylim = c(0,1))
points(1:19, NE.vary.adjR2[[1]], pch = 17, col = "forestgreen")

#NE Aus group 2
plot(1:19, NE.cons.adjR2[[2]], pch = 16, col = "firebrick", ylim = c(0,1))
points(1:19, NE.vary.adjR2[[2]], pch = 17, col = "forestgreen")

#NE Aus group 3
plot(1:19, NE.cons.adjR2[[3]], pch = 16, col = "firebrick", ylim = c(0,1))
points(1:19, NE.vary.adjR2[[3]], pch = 17, col = "forestgreen")

#SE Aus 
SE.cons.adjR2 <- NULL
SE.vary.adjR2 <- NULL
for (j in 1:3) {
  cons.adjR2 <- numeric(19)
  vary.adjR2 <- numeric(19)
  for (i in 1:19) {
    cons.lm.temp <- SErefit.new[[2]][[i]]
    vary.lm.temp <- SErefit.new[[3]][[i]]
    
    cons.adjR2[i] <- summary(cons.lm.temp[[j]])$adj.r.squared
    vary.adjR2[i] <- summary(vary.lm.temp[[j]])$adj.r.squared
  }
  SE.cons.adjR2[[j]] <- cons.adjR2
  SE.vary.adjR2[[j]] <- vary.adjR2
}

#SE Aus group 1
plot(1:19, SE.cons.adjR2[[1]], pch = 16, col = "firebrick", ylim = c(0,1))
points(1:19, SE.vary.adjR2[[1]], pch = 17, col = "forestgreen")

#SE Aus group 2
plot(1:19, SE.cons.adjR2[[2]], pch = 16, col = "firebrick", ylim = c(0,1))
points(1:19, SE.vary.adjR2[[2]], pch = 17, col = "forestgreen")

#SE Aus group 3
plot(1:19, SE.cons.adjR2[[3]], pch = 16, col = "firebrick", ylim = c(0,1))
points(1:19, SE.vary.adjR2[[3]], pch = 17, col = "forestgreen")


## 2019-2020 Season only
#NE Aus
#Group 1
NE.cons.adjR2[[1]][19]
NE.vary.adjR2[[1]][19]
#Group 2
NE.cons.adjR2[[2]][19]
NE.vary.adjR2[[2]][19]
#Group 3
NE.cons.adjR2[[3]][19]
NE.vary.adjR2[[3]][19]

#SE Aus
#Group 1
SE.cons.adjR2[[1]][19]
SE.vary.adjR2[[1]][19]
#Group 2
SE.cons.adjR2[[2]][19]
SE.vary.adjR2[[2]][19]
#Group 3
SE.cons.adjR2[[3]][19]
SE.vary.adjR2[[3]][19]


# . (n.) R^2 (out-of-sample/prediction)
#TODO: update in `model_validation.R`
#TODO: then create some outputs 

# . (n.) K-Fold CV
#TODO: add legend to each plot, then create output (.png)

#NE Aus
NE1.kcv <- lapply(NE.kcv.ebic, function(x) x[[1]])
NE2.kcv <- lapply(NE.kcv.ebic, function(x) x[[2]])
NE3.kcv <- lapply(NE.kcv.ebic, function(x) x[[3]])

NE.kcv.range <- range(NE1.kcv, NE2.kcv, NE3.kcv)
gamma.lab <- sapply(1:12, function(x) paste0("Gamma: ", round(gamma.seq[x], 3)))

plot(1:12, NE1.kcv, type = "b", pch = 16, col = "firebrick3",
     ylab = "K-fold CV", xlab = "",
     ylim = NE.kcv.range, axes = FALSE)
lines(1:12, NE2.kcv, type = "b", pch =17, col = "royalblue4")
lines(1:12, NE3.kcv, type = "b", pch =17, col = "forestgreen")
box()
axis(1, labels = gamma.lab, at = 1:12, cex.axis = 0.67)
axis(2)  
title("NE Aus K-fold CV")

#SE Aus
SE1.kcv <- lapply(SE.kcv.ebic, function(x) x[[1]])
SE2.kcv <- lapply(SE.kcv.ebic, function(x) x[[2]])
SE3.kcv <- lapply(SE.kcv.ebic, function(x) x[[3]])

SE.kcv.range <- range(SE1.kcv, SE2.kcv, SE3.kcv)
gamma.lab <- sapply(1:12, function(x) paste0("Gamma: ", round(gamma.seq[x], 3)))

plot(1:12, SE1.kcv, type = "b", pch = 16, col = "firebrick3",
     ylab = "K-fold CV", xlab = "",
     ylim = SE.kcv.range, axes = FALSE)
lines(1:12, SE2.kcv, type = "b", pch =17, col = "royalblue4")
lines(1:12, SE3.kcv, type = "b", pch =17, col = "forestgreen")
box()
axis(1, labels = gamma.lab, at = 1:12, cex.axis = 0.67)
axis(2)  
title("SE Aus K-fold CV")




# . (number.) RMSE plots
## RMSE setup
NErmse <- NEvalid[[1]]
NErmse.eBIC <- NEvalid.eBIC[[1]]

## -- base : model fits

#NE Aus group 1
NEgroup1.base.rmse <- unlist(lapply(NErmse, function(x) x[[1]][1]))

#NEbase.rmse.eBIC <- matrix(NA, ncol = 19)
k <- 3
NEbase.rmse.eBIC <- unlist(lapply(NErmse.eBIC[[k]], function(x) x[[1]][1]))
k <- 9
NEbase.rmse.eBIC2 <- unlist(lapply(NErmse.eBIC[[k]], function(x) x[[1]][1]))

#boxplot
rmse.range <- range(NEgroup1.base.rmse, NEbase.rmse.eBIC, NEbase.rmse.eBIC2)

model.names <- c("NE Aus 1", paste0("Gamma: ", round(gamma.seq[3], 3)),
                 paste0("Gamma: ", round(gamma.seq[9], 3)))


boxplot(NEgroup1.base.rmse, NEbase.rmse.eBIC, NEbase.rmse.eBIC2,
        pch = 16, 
        ylim = rmse.range, ylab = "RMSE", axes = FALSE)
box()
axis(1, labels = model.names, at = 1:3, cex.axis = 1)
axis(2)                      
title("NE Aus Group 1 Base Model : RMSE", adj = 0)


## --- constant : model fits

#NE Aus Group 1
NEgroup1.const.rmse <- unlist(lapply(NErmse, function(x) x[[2]][1]))
#eBIC setup
k <- 3
NEgroup1.const.rmse.eBIC <- unlist(lapply(NErmse.eBIC[[k]], function(x) x[[2]][1]))
k <- 9
NEgroup1.const.rmse.eBIC2 <- unlist(lapply(NErmse.eBIC[[k]], function(x) x[[2]][1]))

#boxplot
rmse.range <- range(NEgroup1.const.rmse, NEgroup1.const.rmse.eBIC, NEgroup1.const.rmse.eBIC2)

model.names <- c("NE Aus 1", paste0("Gamma: ", round(gamma.seq[3], 3)),
                 paste0("Gamma: ", round(gamma.seq[9], 3)))

#boxplot
boxplot(NEgroup1.const.rmse, NEgroup1.const.rmse.eBIC, NEgroup1.const.rmse.eBIC2,
        pch = 16, 
        ylim = rmse.range, ylab = "RMSE", axes = FALSE)
box()
axis(1, labels = model.names, at = 1:3, cex.axis = 1)
axis(2)                      
title("NE Aus Group 1 Constant Model : RMSE", adj = 0)


#NE Aus Group 2
NEgroup2.const.rmse <- unlist(lapply(NErmse, function(x) x[[2]][2]))
#eBIC setup
k <- 7
NEgroup2.const.rmse.eBIC <- unlist(lapply(NErmse.eBIC[[k]], function(x) x[[2]][2]))

#boxplot
rmse.range <- range(NEgroup2.const.rmse, NEgroup2.const.rmse.eBIC)

model.names <- c("NE Aus 2", paste0("Gamma: ", round(gamma.seq[7], 3)) )

#boxplot
boxplot(NEgroup2.const.rmse, NEgroup2.const.rmse.eBIC,
        pch = 16, 
        ylim = rmse.range, ylab = "RMSE", axes = FALSE)
box()
axis(1, labels = model.names, at = 1:2, cex.axis = 1)
axis(2)                      
title("NE Aus Group 2 Constant Model : RMSE", adj = 0)


#NE Aus Group 3
NEgroup3.const.rmse <- unlist(lapply(NErmse, function(x) x[[2]][3]))

#eBIC setup
k <- 8
NEgroup3.const.rmse.eBIC <- unlist(lapply(NErmse.eBIC[[k]], function(x) x[[2]][3]))

#boxplot
rmse.range <- range(NEgroup3.const.rmse, NEgroup3.const.rmse.eBIC)

model.names <- c("NE Aus 3", paste0("Gamma: ", round(gamma.seq[8], 3)) )


boxplot(NEgroup3.const.rmse, NEgroup3.const.rmse.eBIC,
        pch = 16, 
        ylim = rmse.range, ylab = "RMSE", axes = FALSE)
box()
axis(1, labels = model.names, at = 1:2, cex.axis = 1)
axis(2)                      
title("NE Aus Group 3 Constant Model : RMSE", adj = 0)



# . (number.) CRPS plots
## CPRS setup
NEcprs <- NEvalid[[2]]
NEcprs.eBIC <- NEvalid.eBIC[[2]]


## --- constant : model fits

#NE Aus Group 1
NEgroup1.const.cprs <- unlist(lapply(NEcprs, function(x) x[[2]][1]))
#eBIC setup
k <- 3
NEgroup1.const.cprs.eBIC <- unlist(lapply(NEcprs.eBIC[[k]], function(x) x[[2]][1]))
k <- 9
NEgroup1.const.cprs.eBIC2 <- unlist(lapply(NEcprs.eBIC[[k]], function(x) x[[2]][1]))

#boxplot
cprs.range <- range(NEgroup1.const.cprs, NEgroup1.const.cprs.eBIC, NEgroup1.const.cprs.eBIC2)

model.names <- c("NE Aus 1", paste0("Gamma: ", round(gamma.seq[3], 3)),
                 paste0("Gamma: ", round(gamma.seq[9], 3)))

#boxplot
boxplot(NEgroup1.const.cprs, NEgroup1.const.cprs.eBIC, NEgroup1.const.cprs.eBIC2,
        pch = 16, 
        ylim = cprs.range, ylab = "CRPS", axes = FALSE)
box()
axis(1, labels = model.names, at = 1:3, cex.axis = 1)
axis(2)                      
title("NE Aus Group 1 Constant Model : CRPS", adj = 0)


#NE Aus Group 2
NEgroup2.const.cprs <- unlist(lapply(NEcprs, function(x) x[[2]][2]))
#eBIC setup
k <- 7
NEgroup2.const.cprs.eBIC <- unlist(lapply(NEcprs.eBIC[[k]], function(x) x[[2]][2]))

#boxplot
cprs.range <- range(NEgroup2.const.cprs, NEgroup2.const.cprs.eBIC)

model.names <- c("NE Aus 2", paste0("Gamma: ", round(gamma.seq[7], 3)) )

#boxplot
boxplot(NEgroup2.const.cprs, NEgroup2.const.cprs.eBIC,
        pch = 16, 
        ylim = cprs.range, ylab = "CRPS", axes = FALSE)
box()
axis(1, labels = model.names, at = 1:2, cex.axis = 1)
axis(2)                      
title("NE Aus Group 2 Constant Model : CRPS", adj = 0)


#NE Aus Group 3
NEgroup3.const.cprs <- unlist(lapply(NEcprs, function(x) x[[2]][3]))

#eBIC setup
k <- 8
NEgroup3.const.cprs.eBIC <- unlist(lapply(NEcprs.eBIC[[k]], function(x) x[[2]][3]))

#boxplot
cprs.range <- range(NEgroup3.const.cprs, NEgroup3.const.cprs.eBIC)

model.names <- c("NE Aus 3", paste0("Gamma: ", round(gamma.seq[8], 3)) )


boxplot(NEgroup3.const.cprs, NEgroup3.const.cprs.eBIC,
        pch = 16, 
        ylim = cprs.range, ylab = "CRPS", axes = FALSE)
box()
axis(1, labels = model.names, at = 1:2, cex.axis = 1)
axis(2)                      
title("NE Aus Group 3 Constant Model : CRPS", adj = 0)





# . (number.) Interval Score plots
## ints setup
NEints <- NEvalid[[3]]
NEints.eBIC <- NEvalid.eBIC[[3]]


## --- constant : model fits

#NE Aus Group 1
NEgroup1.const.ints <- unlist(lapply(NEints, function(x) x[[2]][1]))
#eBIC setup
k <- 3
NEgroup1.const.ints.eBIC <- unlist(lapply(NEints.eBIC[[k]], function(x) x[[2]][1]))
k <- 9
NEgroup1.const.ints.eBIC2 <- unlist(lapply(NEints.eBIC[[k]], function(x) x[[2]][1]))

#boxplot
ints.range <- range(NEgroup1.const.ints, NEgroup1.const.ints.eBIC, NEgroup1.const.ints.eBIC2)

model.names <- c("NE Aus 1", paste0("Gamma: ", round(gamma.seq[3], 3)),
                 paste0("Gamma: ", round(gamma.seq[9], 3)))

#boxplot
boxplot(NEgroup1.const.ints, NEgroup1.const.ints.eBIC, NEgroup1.const.ints.eBIC2,
        pch = 16, 
        ylim = ints.range, ylab = "Interval Score", axes = FALSE)
box()
axis(1, labels = model.names, at = 1:3, cex.axis = 1)
axis(2)                      
title("NE Aus Group 1 Constant Model : Interval Score", adj = 0)


#NE Aus Group 2
NEgroup2.const.ints <- unlist(lapply(NEints, function(x) x[[2]][2]))
#eBIC setup
k <- 7
NEgroup2.const.ints.eBIC <- unlist(lapply(NEints.eBIC[[k]], function(x) x[[2]][2]))

#boxplot
ints.range <- range(NEgroup2.const.ints, NEgroup2.const.ints.eBIC)

model.names <- c("NE Aus 2", paste0("Gamma: ", round(gamma.seq[7], 3)) )

#boxplot
boxplot(NEgroup2.const.ints, NEgroup2.const.ints.eBIC,
        pch = 16, 
        ylim = ints.range, ylab = "Interval Score", axes = FALSE)
box()
axis(1, labels = model.names, at = 1:2, cex.axis = 1)
axis(2)                      
title("NE Aus Group 2 Constant Model : Interval Score", adj = 0)


#NE Aus Group 3
NEgroup3.const.ints <- unlist(lapply(NEints, function(x) x[[2]][3]))
#eBIC setup
k <- 8
NEgroup3.const.ints.eBIC <- unlist(lapply(NEints.eBIC[[k]], function(x) x[[2]][3]))

#boxplot
ints.range <- range(NEgroup3.const.ints, NEgroup3.const.ints.eBIC)

model.names <- c("NE Aus 3", paste0("Gamma: ", round(gamma.seq[8], 3)) )


boxplot(NEgroup3.const.ints, NEgroup3.const.ints.eBIC,
        pch = 16, 
        ylim = ints.range, ylab = "Interval Score", axes = FALSE)
box()
axis(1, labels = model.names, at = 1:2, cex.axis = 1)
axis(2)                      
title("NE Aus Group 3 Constant Model : Interval Score", adj = 0)


