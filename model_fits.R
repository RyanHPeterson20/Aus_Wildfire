
#notes:
##-next steps in finalizing the model 
##-add duplicate to finalproject_code folder
##-remove any code and place into refit_cv.rmd when not needed

#TODO: add more libraries as needed
#library 
suppressMessages( library(glmnet)) #test ridge regression for coefs (might not be needed)
suppressMessages( library(RAMP)) #Lasso with efficient solution path.

suppressMessages( library( Metrics)) #measurement metrics

#data import
#TODO: collect data into fewer locations
setwd("~/CO_AUS/Aus_CO-main/Interactions")
load( "Data/ne_data.rda")
load( "Data/se_data.rda")
load( "Data/bounded_data.rda")
load( "Data/data_matrix.rda")
load( "Data/lag_list.rda")
load( "Data/data_quantile.rda")

setwd("~/CO_AUS/Aus_CO-main/Interactions_New")
load("IOD_lag.rda") #IOD data

#function load
#TODO: organize and combine helper functions
setwd("~/CO_AUS/Aus_CO-main/Interactions")
#source("group_functions.R") #grouping/clustering
source("group_functionsNew.R") #new groupings
source("refit_functions.R") #coef/refit functions (for hiernet, currently not used)

setwd("~/CO_AUS/Aus_CO-main/Interactions_New")
source("predictionTest_functions.R") #predtest functions 
source("refit_helper.R") #RAMP refit functions

# Setup
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

#full model
NEpreds_new <- NElag_3group(NE_laglist = NE_iodlag, j = 1:19)
NEresp_new <- NEresp_3group(NEAus_mat = NEAus_mat, j = 1:19)

SEpreds_new <- SElag_3group(SE_laglist = SE_iodlag, j = 1:19)
SEresp_new <- SEresp_3group(SEAus_mat = SEAus_mat, j = 1:19)

##-----Main-----##

##NE Aus Models with RAMP
#NE Aus Group 1
y.1 <- as.numeric(NEresp_new[[1]]) #co response
X.1 <- cbind(as.matrix(NEpreds_new[[1]][ ,1:312])) #preds with OLR

NE1.ramp <- RAMP(X = X.1, y = y.1,
                 penalty = "LASSO",
                 tune = "BIC",
                 n.lambda = 500)

#terms only
NE1.terms <- terms_only(NE1.ramp, X.1)
NE1.refit <- refit_ramp(NE1.ramp, X.1)
#lm
lm.data.1 <- as.data.frame(cbind(y.1, X.1))
names(lm.data.1)[1] <- "co"

NE1.lm <- lm(formula(NE1.refit), lm.data.1)

#NE Aus Group 2
y.2 <- as.numeric(NEresp_new[[2]]) #co response
X.2 <- cbind(as.matrix(NEpreds_new[[2]][ ,1:312])) #preds with OLR

NE2.ramp <- RAMP(X = X.2, y = y.2,
                 penalty = "LASSO",
                 tune = "BIC",
                 n.lambda = 500)

#terms only
NE2.terms <- terms_only(NE2.ramp, X.2)
NE2.refit <- refit_ramp(NE2.ramp, X.2)
#lm
lm.data.2 <- as.data.frame(cbind(y.2, X.2))
names(lm.data.2)[1] <- "co"

NE2.lm <- lm(formula(NE2.refit), lm.data.2)

#NE Aus Group 3
y.3 <- as.numeric(NEresp_new[[3]]) #co response
X.3 <- cbind(as.matrix(NEpreds_new[[3]][ ,1:312])) #preds with OLR

NE3.ramp <- RAMP(X = X.3, y = y.3,
                 penalty = "LASSO",
                 tune = "BIC",
                 n.lambda = 500)

#terms only
NE3.terms <- terms_only(NE3.ramp, X.3)
NE3.refit <- refit_ramp(NE3.ramp, X.3)
#lm
lm.data.3 <- as.data.frame(cbind(y.3, X.3))
names(lm.data.3)[1] <- "co"

NE3.lm <- lm(formula(NE3.refit), lm.data.3)


##SE Aus Models with RAMP
#SE Aus Group 1
y.1 <- as.numeric(SEresp_new[[1]]) #co response
X.1 <- cbind(as.matrix(SEpreds_new[[1]][ ,1:312])) #preds with OLR

SE1.ramp <- RAMP(X = X.1, y = y.1,
                 penalty = "LASSO",
                 tune = "BIC",
                 n.lambda = 500)

#terms only
SE1.terms <- terms_only(SE1.ramp, X.1)
SE1.refit <- refit_ramp(SE1.ramp, X.1)
#lm
lm.data.1 <- as.data.frame(cbind(y.1, X.1))
names(lm.data.1)[1] <- "co"

SE1.lm <- lm(formula(SE1.refit), lm.data.1)

#SE Aus Group 2
y.2 <- as.numeric(SEresp_new[[2]]) #co response
X.2 <- cbind(as.matrix(SEpreds_new[[2]][ ,1:312])) #preds with OLR

SE2.ramp <- RAMP(X = X.2, y = y.2,
                 penalty = "LASSO",
                 tune = "BIC",
                 n.lambda = 500)

#terms only
SE2.terms <- terms_only(SE2.ramp, X.2)
SE2.refit <- refit_ramp(SE2.ramp, X.2)
#lm
lm.data.2 <- as.data.frame(cbind(y.2, X.2))
names(lm.data.2)[1] <- "co"

SE2.lm <- lm(formula(SE2.refit), lm.data.2)

#SE Aus Group 3
y.3 <- as.numeric(SEresp_new[[3]]) #co response
X.3 <- cbind(as.matrix(SEpreds_new[[3]][ ,1:312])) #preds with OLR

SE3.ramp <- RAMP(X = X.3, y = y.3,
                 penalty = "LASSO",
                 tune = "BIC",
                 n.lambda = 500)

#terms only
SE3.terms <- terms_only(SE3.ramp, X.3)
SE3.refit <- refit_ramp(SE3.ramp, X.3)
#lm
lm.data.3 <- as.data.frame(cbind(y.3, X.3))
names(lm.data.3)[1] <- "co"

SE3.lm <- lm(formula(SE3.refit), lm.data.3)

#base model fit output
NEmodels <- list(NE1.lm, NE2.lm, NE3.lm)
NErefits <- list(NE1.refit, NE2.refit, NE3.refit)
SEmodels <- list(SE1.lm, SE2.lm, SE3.lm)
SErefits <- list(SE1.refit, SE2.refit, SE3.refit)

setwd("~/CO_AUS/Aus_CO-main/Interactions_New")
save(NEmodels, SEmodels, 
     NErefits, SErefits, file = "base_RAMPmodels.rda")


##eBIC Model Fits
#gamma
gamma.seq <- seq(0, 1, length.out = 12)

#NE Aus Group 1
y.1 <- as.numeric(NEresp_new[[1]]) #co response
X.1 <- cbind(as.matrix(NEpreds_new[[1]][ ,1:312])) #preds with OLR

NE1.eBIC.lm <- list()
NE1.eBIC.refit <- list()
for (j in 1:length(gamma.seq)) {
  negroup1.temp <- RAMP(X = X.1, y = y.1,
                        penalty = "LASSO",
                        tune = "EBIC",
                        ebic.gamma = gamma.seq[j],
                        n.lambda = 500)
  
  #terms only
  NE1.eBIC.temp <- refit_ramp(negroup1.temp, X.1)
  NE1.eBIC.refit[[paste0("Gamma: ", round(gamma.seq[[j]], 3))]] <- NE1.eBIC.temp
  
  #refit with lm 
  lm.data.1 <- as.data.frame(cbind(y.1, X.1))
  names(lm.data.1)[1] <- "co"
  
  NE1.eBIC.lm[[paste0("Gamma: ", round(gamma.seq[[j]], 3))]] <- lm(formula(NE1.eBIC.temp), lm.data.1)
}


#NE Aus Group 2
y.2 <- as.numeric(NEresp_new[[2]]) #co response
X.2 <- cbind(as.matrix(NEpreds_new[[2]][ ,1:312])) #preds with OLR

NE2.eBIC.lm <- list()
NE2.eBIC.refit <- list()
for (j in 1:length(gamma.seq)) {
  negroup2.temp <- RAMP(X = X.2, y = y.2,
                        penalty = "LASSO",
                        tune = "EBIC",
                        ebic.gamma = gamma.seq[j],
                        n.lambda = 500)
  
  #terms only
  NE2.eBIC.temp <- refit_ramp(negroup2.temp, X.2)
  NE2.eBIC.refit[[paste0("Gamma: ", round(gamma.seq[[j]], 3))]] <- NE2.eBIC.temp
  
  #refit with lm 
  lm.data.2 <- as.data.frame(cbind(y.2, X.2))
  names(lm.data.2)[1] <- "co"
  
  NE2.eBIC.lm[[paste0("Gamma: ", round(gamma.seq[[j]], 3))]] <- lm(formula(NE2.eBIC.temp), lm.data.2)
}


#NE Aus Group 3
y.3 <- as.numeric(NEresp_new[[3]]) #co response
X.3 <- cbind(as.matrix(NEpreds_new[[3]][ ,1:312])) #preds with OLR

NE3.eBIC.lm <- list()
NE3.eBIC.refit <- list()
for (j in 1:length(gamma.seq)) {
  negroup3.temp <- RAMP(X = X.3, y = y.3,
                        penalty = "LASSO",
                        tune = "EBIC",
                        ebic.gamma = gamma.seq[j],
                        n.lambda = 500)
  
  #terms only
  NE3.eBIC.temp <- refit_ramp(negroup3.temp, X.3)
  NE3.eBIC.refit[[paste0("Gamma: ", round(gamma.seq[[j]], 3))]] <- NE3.eBIC.temp
  
  #refit with lm 
  lm.data.3 <- as.data.frame(cbind(y.3, X.3))
  names(lm.data.3)[1] <- "co"
  
  NE3.eBIC.lm[[paste0("Gamma: ", round(gamma.seq[[j]], 3))]] <- lm(formula(NE3.eBIC.temp), lm.data.3)
}



#SE Aus Group 1
y.1 <- as.numeric(SEresp_new[[1]]) #co response
X.1 <- cbind(as.matrix(SEpreds_new[[1]][ ,1:312])) #preds with OLR

SE1.eBIC.lm <- list()
SE1.eBIC.refit <- list()
for (j in 1:length(gamma.seq)) {
  segroup1.temp <- RAMP(X = X.1, y = y.1,
                        penalty = "LASSO",
                        tune = "EBIC",
                        ebic.gamma = gamma.seq[j],
                        n.lambda = 500)
  
  #terms only
  SE1.eBIC.temp <- refit_ramp(segroup1.temp, X.1)
  SE1.eBIC.refit[[paste0("Gamma: ", round(gamma.seq[[j]], 3))]] <- SE1.eBIC.temp
  
  #refit with lm 
  lm.data.1 <- as.data.frame(cbind(y.1, X.1))
  names(lm.data.1)[1] <- "co"
  
  SE1.eBIC.lm[[paste0("Gamma: ", round(gamma.seq[[j]], 3))]] <- lm(formula(SE1.eBIC.temp), lm.data.1)
}


#SE Aus Group 2
y.2 <- as.numeric(SEresp_new[[2]]) #co response
X.2 <- cbind(as.matrix(SEpreds_new[[2]][ ,1:312])) #preds with OLR

SE2.eBIC.lm <- list()
SE2.eBIC.refit <- list()
for (j in 1:length(gamma.seq)) {
  segroup2.temp <- RAMP(X = X.2, y = y.2,
                        penalty = "LASSO",
                        tune = "EBIC",
                        ebic.gamma = gamma.seq[j],
                        n.lambda = 500)
  
  #terms only
  SE2.eBIC.temp <- refit_ramp(segroup2.temp, X.2)
  SE2.eBIC.refit[[paste0("Gamma: ", round(gamma.seq[[j]], 3))]] <- SE2.eBIC.temp
  
  #refit with lm 
  lm.data.2 <- as.data.frame(cbind(y.2, X.2))
  names(lm.data.2)[1] <- "co"
  
  SE2.eBIC.lm[[paste0("Gamma: ", round(gamma.seq[[j]], 3))]] <- lm(formula(SE2.eBIC.temp), lm.data.2)
}


#SE Aus Group 3
y.3 <- as.numeric(SEresp_new[[3]]) #co response
X.3 <- cbind(as.matrix(SEpreds_new[[3]][ ,1:312])) #preds with OLR

SE3.eBIC.lm <- list()
SE3.eBIC.refit <- list()
for (j in 1:length(gamma.seq)) {
  segroup3.temp <- RAMP(X = X.3, y = y.3,
                        penalty = "LASSO",
                        tune = "EBIC",
                        ebic.gamma = gamma.seq[j],
                        n.lambda = 500)
  
  #terms only
  SE3.eBIC.temp <- refit_ramp(segroup3.temp, X.3)
  SE3.eBIC.refit[[paste0("Gamma: ", round(gamma.seq[[j]], 3))]] <- SE3.eBIC.temp
  
  #refit with lm 
  lm.data.3 <- as.data.frame(cbind(y.3, X.3))
  names(lm.data.3)[1] <- "co"
  
  SE3.eBIC.lm[[paste0("Gamma: ", round(gamma.seq[[j]], 3))]] <- lm(formula(SE3.eBIC.temp), lm.data.3)
}


#output eBIC models
NEmodels.eBIC <- list(NE1.eBIC.lm, NE2.eBIC.lm, NE3.eBIC.lm)
NErefits.eBIC <- list(NE1.eBIC.refit, NE2.eBIC.refit, NE3.eBIC.refit)
SEmodels.eBIC <- list(SE1.eBIC.lm, SE2.eBIC.lm, SE3.eBIC.lm)
SErefits.eBIC <- list(SE1.eBIC.refit, SE2.eBIC.refit, SE3.eBIC.refit)

setwd("~/CO_AUS/Aus_CO-main/Interactions_New")
save(NEmodels.eBIC, SEmodels.eBIC, 
     NErefits.eBIC, SErefits.eBIC, file = "eBIC_RAMPmodels.rda")

