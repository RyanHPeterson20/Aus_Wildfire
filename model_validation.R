#model validation

##-include model refits (e.g. LOO)
###--include RMSE, CPRS, and Int.Score
##-perform refit and new RAMP for k-fold CV

#TODO: import saved fits from model_fits.R
##-format LOO (leave-one year-out) loops

# libraries
#suppressMessages( library())
suppressMessages( library(RAMP)) #Lasso with efficient solution path.

suppressMessages( library( Metrics)) #measurement metrics

#parallelization setup
suppressMessages( library(foreach)) 
suppressMessages( library(parallel))
suppressMessages( library(doParallel))

# data load
#TODO: update base data into a fewer files
## and combine locations (this will be done in final project in some abstract Data folder)
setwd("~/CO_AUS/Aus_CO-main/Interactions")
load( "Data/ne_data.rda")
load( "Data/se_data.rda")
load( "Data/bounded_data.rda")
load( "Data/data_matrix.rda")
load( "Data/lag_list.rda")

setwd("~/CO_AUS/Aus_CO-main/Interactions_New")
load("IOD_lag.rda") #IOD data

#model fits from `model_fits.R`
setwd("~/CO_AUS/Aus_CO-main/Interactions_New")
load("base_RAMPmodels.rda") # 'base' fits (BIC)
load("eBIC_RAMPmodels.rda") # fits with eBIC

# function load
#TODO: organize and combine helper functions
setwd("~/CO_AUS/Aus_CO-main/Interactions")
#source("group_functions.R") #grouping/clustering
source("group_functionsNew.R") #new groupings
source("refit_functions.R") #coef/refit functions (for hiernet, currently not used)

setwd("~/CO_AUS/Aus_CO-main/Interactions_New")
source("predictionTest_functions.R") #prediction test functions 
source("refit_helper.R") #RAMP refit functions


## --- setup --- ##

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

#partial model
NEpreds_newpart <- NElag_3group(NE_laglist = NE_iodlag, j = -c(19))
NEresp_newpart <- NEresp_3group(NEAus_mat = NEAus_mat, j = -c(19))

SEpreds_newpart <- SElag_3group(SE_laglist = SE_iodlag, j = -c(19))
SEresp_newpart <- SEresp_3group(SEAus_mat = SEAus_mat, j = -c(19))

## --- partial data setup 
#leaving out a single year (loo)
NEpartial_preds <- list()
NEpartial_resp <- list()
SEpartial_preds <- list()
SEpartial_resp <- list()
for (j in 1:length(seasons)) {
  #NE Aus
  NEpartial_preds[[seasons[j]]] <- NElag_3group(NE_laglist = NE_iodlag, j = -c(j))
  NEpartial_resp[[seasons[j]]] <- NEresp_3group(NEAus_mat = NEAus_mat, j = -c(j))
  
  #SE Aus
  SEpartial_preds[[seasons[j]]] <- SElag_3group(SE_laglist = SE_iodlag, j = -c(j))
  SEpartial_resp[[seasons[j]]] <- SEresp_3group(SEAus_mat = SEAus_mat, j = -c(j))
}

#extracting a single year (validation data)
NEvalid_preds <- list()
NEvalid_resp <- list()
SEvalid_preds <- list()
SEvalid_resp <- list()
for (k in 1:length(seasons)) {
  #NE Aus
  NEvalid_preds[[seasons[k]]] <- NElag_3group(NE_laglist = NE_iodlag, j = c(k))
  NEvalid_resp[[seasons[k]]] <- NEresp_3group(NEAus_mat = NEAus_mat, j = c(k))
  
  #SE Aus
  SEvalid_preds[[seasons[k]]] <- SElag_3group(SE_laglist = SE_iodlag, j = c(k))
  SEvalid_resp[[seasons[k]]] <- SEresp_3group(SEAus_mat = SEAus_mat, j = c(k))
}


## --- Main --- ##
#model refits and validation

## - Leave One year Out (LOO)

##NE Aus refits and validation
#model validation
NE.rmse <- NULL
NE.cprs <- NULL
NE.ints <- NULL
NE.pred.int <- NULL
#models (refit and lm)
NE.vary.terms <- NULL
NE.const.LM <- NULL 
NE.vary.LM <- NULL
for (i in 1:length(seasons)) {
  #data w/o season (train)
  fit.resp <- NEpartial_resp[[i]]
  fit.preds <- NEpartial_preds[[i]]
  
  #data w/ only season (test)
  valid.resp <- NEvalid_resp[[i]]
  valid.preds <- NEvalid_preds[[i]]  
  
  #group data objects
  NE.var.refit <- NULL #varying terms
  NE.con <- NULL #constant linear models
  NE.var <- NULL #varying linear models
  NErmse.yearly <- matrix(NA, ncol = 3)
  colnames(NErmse.yearly) <- c("base.pred", "const.pred", "vary.pred")
  NEcprs.yearly <- matrix(NA, ncol = 3)
  colnames(NEcprs.yearly) <- c("base.pred", "const.pred", "vary.pred")
  NEis.yearly <- matrix(NA, ncol = 3)
  colnames(NEis.yearly) <- c("base.pred", "const.pred", "vary.pred")
  NE.intervals <- matrix(NA, ncol = 10)
  colnames(NE.intervals) <- c("true", "base.fit", "base.lwr", "base.upr",  
                              "const.fit", "const.lwr", "const.upr",
                              "vary.fit", "vary.lwr", "vary.upr")
  for (j in 1:3) {
    #get base model terms (and fits)
    NE.lm.base <- NEmodels[[j]] #lm model for NE group j 
    NE.terms.base <- NErefits[[j]] #terms for NE group j
    
    #lm fit data setup
    y.fit <- as.numeric(fit.resp[[j]])
    #with OLR
    X.fit <- cbind(as.matrix(fit.preds[[j]][ ,1:312]))
    
    lm.data.fit <- as.data.frame(cbind(y.fit, X.fit))
    names(lm.data.fit)[1] <- "co"
    
    #varying ramp fit
    vary.fit <- RAMP(X = X.fit, y = y.fit,
                     penalty = "LASSO",
                     tune = "BIC",
                     n.lambda = 500)
    
    NE.terms.vary <- refit_ramp(vary.fit, X.fit)
    
    #refit
    NE.lm.const <- lm(formula(NE.terms.base), lm.data.fit)
    NE.lm.vary <- lm(formula(NE.terms.vary), lm.data.fit)
    
    #assign terms and models
    NE.var.refit[[j]] <- NE.terms.vary
    NE.con[[j]] <- NE.lm.const
    NE.var[[j]] <- NE.lm.vary
    
    #prediction and validation 
    y.valid <- as.numeric(valid.resp[[j]])
    X.valid <- valid.preds[[j]][ ,1:312]
    
    pred.base <- predict(NE.lm.base, X.valid, se.fit = TRUE)
    pred.const <- predict(NE.lm.const, X.valid, se.fit = TRUE)
    pred.vary <- predict(NE.lm.vary, X.valid, se.fit = TRUE)
    
    #rmse
    rmse.base <-  rmse(y.valid, pred.base$fit)
    rmse.const <- rmse(y.valid, pred.const$fit)
    rmse.vary <- rmse(y.valid, pred.vary$fit)
    #cprs
    cprs.base <- mean(CPRS(list(mean = pred.base$fit, sd = pred.base$se.fit), y.valid))
    cprs.const <- mean(CPRS(list(mean = pred.const$fit, sd = pred.const$se.fit), y.valid))
    cprs.vary <- mean(CPRS(list(mean = pred.vary$fit, sd = pred.vary$se.fit), y.valid))
    #IS95
    is.base <- mean(intscore(list(mean = pred.base$fit, sd = pred.base$se.fit), y.valid))
    is.const <- mean(intscore(list(mean = pred.const$fit, sd = pred.const$se.fit), y.valid))
    is.vary <- mean(intscore(list(mean = pred.vary$fit, sd = pred.vary$se.fit), y.valid))    
    
    #assign validations
    NErmse.yearly <- rbind(NErmse.yearly, cbind(rmse.base, rmse.const, rmse.vary))
    NEcprs.yearly <- rbind(NEcprs.yearly, cbind(cprs.base, cprs.const, cprs.vary))
    NEis.yearly <- rbind(NEis.yearly, cbind(is.base, is.const, is.vary))
    
    #intervals
    pred.base.interval <- predict(NE.lm.base, X.valid, interval = "prediction")
    pred.const.interval <- predict(NE.lm.const, X.valid, interval = "prediction")
    pred.vary.interval <- predict(NE.lm.vary, X.valid, interval = "prediction")
    
    #assign intervals
    NE.intervals <- rbind(NE.intervals, 
                          cbind(y.valid, pred.base.interval, pred.const.interval, pred.vary.interval))
  }
  NE.vary.terms[[seasons[i]]] <- NE.var.refit
  NE.const.LM[[seasons[i]]] <- NE.con
  NE.vary.LM[[seasons[i]]] <- NE.var
  
  NE.rmse[[seasons[i]]] <- as.data.frame(NErmse.yearly[-1, ])
  NE.cprs[[seasons[i]]] <- as.data.frame(NEcprs.yearly[-1, ])
  NE.ints[[seasons[i]]] <- as.data.frame(NEis.yearly[-1, ])
  NE.pred.int[[seasons[i]]] <- as.data.frame(NE.intervals[-1, ])
  
}

NEvalid <- list(rmse = NE.rmse, cprs = NE.cprs, ints = NE.ints, interval = NE.pred.int)
NErefit.new <- list(NE.vary.terms, NE.const.LM, NE.vary.LM)

##SE Aus refits and validation
#model validation
SE.rmse <- NULL
SE.cprs <- NULL
SE.ints <- NULL
SE.pred.int <- NULL
#models (refit and lm)
SE.vary.terms <- NULL
SE.const.LM <- NULL 
SE.vary.LM <- NULL
for (i in 1:length(seasons)) {
  #data w/o season (train)
  fit.resp <- SEpartial_resp[[i]]
  fit.preds <- SEpartial_preds[[i]]
  
  #data w/ only season (test)
  valid.resp <- SEvalid_resp[[i]]
  valid.preds <- SEvalid_preds[[i]]  
  
  #group data objects
  SE.var.refit <- NULL #varying terms
  SE.con <- NULL #constant linear models
  SE.var <- NULL #varying linear models
  SErmse.yearly <- matrix(NA, ncol = 3)
  colnames(SErmse.yearly) <- c("base.pred", "const.pred", "vary.pred")
  SEcprs.yearly <- matrix(NA, ncol = 3)
  colnames(SEcprs.yearly) <- c("base.pred", "const.pred", "vary.pred")
  SEis.yearly <- matrix(NA, ncol = 3)
  colnames(SEis.yearly) <- c("base.pred", "const.pred", "vary.pred")
  SE.intervals <- matrix(NA, ncol = 10)
  colnames(SE.intervals) <- c("true", "base.fit", "base.lwr", "base.upr",  
                              "const.fit", "const.lwr", "const.upr",
                              "vary.fit", "vary.lwr", "vary.upr")
  for (j in 1:3) {
    #get base model terms (and fits)
    SE.lm.base <- SEmodels[[j]] #lm model for SE group j 
    SE.terms.base <- SErefits[[j]] #terms for SE group j
    
    #lm fit data setup
    y.fit <- as.numeric(fit.resp[[j]])
    #with OLR
    X.fit <- cbind(as.matrix(fit.preds[[j]][ ,1:312]))
    
    lm.data.fit <- as.data.frame(cbind(y.fit, X.fit))
    names(lm.data.fit)[1] <- "co"
    
    #varying ramp fit
    vary.fit <- RAMP(X = X.fit, y = y.fit,
                     penalty = "LASSO",
                     tune = "BIC",
                     n.lambda = 500)
    
    SE.terms.vary <- refit_ramp(vary.fit, X.fit)
    
    #refit
    SE.lm.const <- lm(formula(SE.terms.base), lm.data.fit)
    SE.lm.vary <- lm(formula(SE.terms.vary), lm.data.fit)
    
    #assign terms and models
    SE.var.refit[[j]] <- SE.terms.vary
    SE.con[[j]] <- SE.lm.const
    SE.var[[j]] <- SE.lm.vary
    
    #prediction and validation 
    y.valid <- as.numeric(valid.resp[[j]])
    X.valid <- valid.preds[[j]][ ,1:312]
    
    pred.base <- predict(SE.lm.base, X.valid, se.fit = TRUE)
    pred.const <- predict(SE.lm.const, X.valid, se.fit = TRUE)
    pred.vary <- predict(SE.lm.vary, X.valid, se.fit = TRUE)
    
    #rmse
    rmse.base <-  rmse(y.valid, pred.base$fit)
    rmse.const <- rmse(y.valid, pred.const$fit)
    rmse.vary <- rmse(y.valid, pred.vary$fit)
    #cprs
    cprs.base <- mean(CPRS(list(mean = pred.base$fit, sd = pred.base$se.fit), y.valid))
    cprs.const <- mean(CPRS(list(mean = pred.const$fit, sd = pred.const$se.fit), y.valid))
    cprs.vary <- mean(CPRS(list(mean = pred.vary$fit, sd = pred.vary$se.fit), y.valid))
    #IS95
    is.base <- mean(intscore(list(mean = pred.base$fit, sd = pred.base$se.fit), y.valid))
    is.const <- mean(intscore(list(mean = pred.const$fit, sd = pred.const$se.fit), y.valid))
    is.vary <- mean(intscore(list(mean = pred.vary$fit, sd = pred.vary$se.fit), y.valid))    
    
    #assign validations
    SErmse.yearly <- rbind(SErmse.yearly, cbind(rmse.base, rmse.const, rmse.vary))
    SEcprs.yearly <- rbind(SEcprs.yearly, cbind(cprs.base, cprs.const, cprs.vary))
    SEis.yearly <- rbind(SEis.yearly, cbind(is.base, is.const, is.vary))
    
    #intervals
    pred.base.interval <- predict(SE.lm.base, X.valid, interval = "prediction")
    pred.const.interval <- predict(SE.lm.const, X.valid, interval = "prediction")
    pred.vary.interval <- predict(SE.lm.vary, X.valid, interval = "prediction")
    
    #assign intervals
    SE.intervals <- rbind(SE.intervals, 
                          cbind(y.valid, pred.base.interval, pred.const.interval, pred.vary.interval))
  }
  SE.vary.terms[[seasons[i]]] <- SE.var.refit
  SE.const.LM[[seasons[i]]] <- SE.con
  SE.vary.LM[[seasons[i]]] <- SE.var
  
  SE.rmse[[seasons[i]]] <- as.data.frame(SErmse.yearly[-1, ])
  SE.cprs[[seasons[i]]] <- as.data.frame(SEcprs.yearly[-1, ])
  SE.ints[[seasons[i]]] <- as.data.frame(SEis.yearly[-1, ])
  SE.pred.int[[seasons[i]]] <- as.data.frame(SE.intervals[-1, ])
}


SEvalid <- list(SE.rmse, SE.cprs, SE.ints, SE.pred.int)
SErefit.new <- list(SE.vary.terms, SE.const.LM, SE.vary.LM)

setwd("~/CO_AUS/Aus_CO-main/Interactions_New")
save(NEvalid, NErefit.new, 
     SEvalid, SErefit.new, file = "validation_refits.rda")


#eBIC refits
##parallel interior loop (model refit for eBIC varying)

#gamma (for eBIC)
gamma.seq <- seq(0, 1, length.out = 12)

#fit first
n.cores <- 10
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)

NE.eBIC.vary <- NULL 
for (i in 1:length(seasons)) { #
  #data w/o season (train)
  fit.resp <- NEpartial_resp[[i]]
  fit.preds <- NEpartial_preds[[i]]

  eBIC.vary.seasonal <- NULL
  for (j in 1:3) {
  
    #lm fit data setup
    y.fit <- as.numeric(fit.resp[[j]])
    #with OLR
    X.fit <- cbind(as.matrix(fit.preds[[j]][ ,1:312]))
    
    #for each model fit
    eBIC.vary.fit <- foreach(k = 1:length(gamma.seq)) %dopar% {
      
      #varying ramp fit (eBIC)
      RAMP::RAMP(X = X.fit, y = y.fit,
                 penalty = "LASSO",
                 tune = "EBIC",
                 ebic.gamma = gamma.seq[k],
                 n.lambda = 500)
    }
    eBIC.vary.seasonal[[j]] <- eBIC.vary.fit

  }
  NE.eBIC.vary[[seasons[i]]] <- eBIC.vary.seasonal
  
}

parallel::stopCluster(cl = my.cluster)  


setwd("~/CO_AUS/Aus_CO-main/Interactions_New")
save(NE.eBIC.vary, file = "tempNEeBIC.rda")

load("tempNEeBIC.rda")

#NE Aus eBIC refits
NE.rmse.eBIC <- NULL
NE.cprs.eBIC <- NULL
NE.ints.eBIC <- NULL
NE.pred.int.eBIC <- NULL
#models (refit and lm)
NE.vary.terms.eBIC <- NULL
NE.const.LM.eBIC <- NULL 
NE.vary.LM.eBIC <- NULL
#set outer loop for gamma values
for (k in 1:length(gamma.seq)) {
  #internal setup
  #model validation
  NE.rmse <- NULL
  NE.cprs <- NULL
  NE.ints <- NULL
  NE.pred.int <- NULL
  #models (refit and lm)
  NE.vary.terms <- NULL
  NE.const.LM <- NULL 
  NE.vary.LM <- NULL
  for (i in 1:length(seasons)) {
    #data w/o season (train)
    fit.resp <- NEpartial_resp[[i]]
    fit.preds <- NEpartial_preds[[i]]
    
    #data w/ only season (test)
    valid.resp <- NEvalid_resp[[i]]
    valid.preds <- NEvalid_preds[[i]]  
    
    #group data objects
    NE.var.refit <- NULL #varying terms
    NE.con <- NULL #constant linear models
    NE.var <- NULL #varying linear models
    NErmse.yearly <- matrix(NA, ncol = 3)
    colnames(NErmse.yearly) <- c("base.pred", "const.pred", "vary.pred")
    NEcprs.yearly <- matrix(NA, ncol = 3)
    colnames(NEcprs.yearly) <- c("base.pred", "const.pred", "vary.pred")
    NEis.yearly <- matrix(NA, ncol = 3)
    colnames(NEis.yearly) <- c("base.pred", "const.pred", "vary.pred")
    NE.intervals <- matrix(NA, ncol = 10)
    colnames(NE.intervals) <- c("true", "base.fit", "base.lwr", "base.upr",  
                                "const.fit", "const.lwr", "const.upr",
                                "vary.fit", "vary.lwr", "vary.upr")
    for (j in 1:3) {
      #get base model terms (and fits)
      NE.lm.base <- NEmodels.eBIC[[j]][[k]] #lm model for NE group j 
      NE.terms.base <- NErefits.eBIC[[j]][[k]] #terms for NE group j
      NE.vary.fit <- NE.eBIC.vary[[i]][[j]][[k]]
      
      #lm fit data setup
      y.fit <- as.numeric(fit.resp[[j]])
      #with OLR
      X.fit <- cbind(as.matrix(fit.preds[[j]][ ,1:312]))
      
      lm.data.fit <- as.data.frame(cbind(y.fit, X.fit))
      names(lm.data.fit)[1] <- "co"
      
      
      NE.terms.vary <- refit_ramp(NE.vary.fit, X.fit)
      
      #refit
      NE.lm.const <- lm(formula(NE.terms.base), lm.data.fit)
      NE.lm.vary <- lm(formula(NE.terms.vary), lm.data.fit)
      
      #assign terms and models
      NE.var.refit[[j]] <- NE.terms.vary
      NE.con[[j]] <- NE.lm.const
      NE.var[[j]] <- NE.lm.vary
      
      #prediction and validation 
      y.valid <- as.numeric(valid.resp[[j]])
      X.valid <- valid.preds[[j]][ ,1:312]
      
      pred.base <- predict(NE.lm.base, X.valid, se.fit = TRUE)
      pred.const <- predict(NE.lm.const, X.valid, se.fit = TRUE)
      pred.vary <- predict(NE.lm.vary, X.valid, se.fit = TRUE)
      
      #rmse
      rmse.base <-  rmse(y.valid, pred.base$fit)
      rmse.const <- rmse(y.valid, pred.const$fit)
      rmse.vary <- rmse(y.valid, pred.vary$fit)
      #cprs
      cprs.base <- mean(CPRS(list(mean = pred.base$fit, sd = pred.base$se.fit), y.valid))
      cprs.const <- mean(CPRS(list(mean = pred.const$fit, sd = pred.const$se.fit), y.valid))
      cprs.vary <- mean(CPRS(list(mean = pred.vary$fit, sd = pred.vary$se.fit), y.valid))
      #IS95
      is.base <- mean(intscore(list(mean = pred.base$fit, sd = pred.base$se.fit), y.valid))
      is.const <- mean(intscore(list(mean = pred.const$fit, sd = pred.const$se.fit), y.valid))
      is.vary <- mean(intscore(list(mean = pred.vary$fit, sd = pred.vary$se.fit), y.valid))    
      
      #assign validations
      NErmse.yearly <- rbind(NErmse.yearly, cbind(rmse.base, rmse.const, rmse.vary))
      NEcprs.yearly <- rbind(NEcprs.yearly, cbind(cprs.base, cprs.const, cprs.vary))
      NEis.yearly <- rbind(NEis.yearly, cbind(is.base, is.const, is.vary))
      
      #intervals
      pred.base.interval <- predict(NE.lm.base, X.valid, interval = "prediction")
      pred.const.interval <- predict(NE.lm.const, X.valid, interval = "prediction")
      pred.vary.interval <- predict(NE.lm.vary, X.valid, interval = "prediction")
      
      #assign intervals
      NE.intervals <- rbind(NE.intervals, 
                            cbind(y.valid, pred.base.interval, pred.const.interval, pred.vary.interval))
    }
    NE.vary.terms[[seasons[i]]] <- NE.var.refit
    NE.const.LM[[seasons[i]]] <- NE.con
    NE.vary.LM[[seasons[i]]] <- NE.var
    
    NE.rmse[[seasons[i]]] <- as.data.frame(NErmse.yearly[-1, ])
    NE.cprs[[seasons[i]]] <- as.data.frame(NEcprs.yearly[-1, ])
    NE.ints[[seasons[i]]] <- as.data.frame(NEis.yearly[-1, ])
    NE.pred.int[[seasons[i]]] <- as.data.frame(NE.intervals[-1, ])
    
  }
  #attach the internals to eBIC lists
  NE.vary.terms.eBIC[[k]] <- NE.vary.terms
  NE.const.LM.eBIC[[k]] <- NE.const.LM 
  NE.vary.LM.eBIC[[k]]  <- NE.vary.LM
  
  NE.rmse.eBIC[[k]]  <- NE.rmse
  NE.cprs.eBIC[[k]]  <- NE.cprs
  NE.ints.eBIC[[k]]  <- NE.ints
  NE.pred.int.eBIC[[k]]  <- NE.pred.int
}

NEvalid.eBIC <- list(NE.rmse.eBIC, NE.cprs.eBIC, NE.ints.eBIC, NE.pred.int.eBIC)
NErefit.new.eBIC <- list(NE.vary.terms.eBIC, NE.const.LM.eBIC, NE.vary.LM.eBIC)

#TODO: repeat above (varying fits and refits) but with SE Aus 



setwd("~/CO_AUS/Aus_CO-main/Interactions_New")
save(NEvalid.eBIC, NErefit.new.eBIC, 
     file = "validation_refitsEBIC.rda")


#perform k-fold cv on both BIC and eBIC
k <- 10

#NE aus
NE.kcv <- NULL
for(j in 1:3){ 
  temp.pred <- NEpreds_new[[j]]
  temp.resp <- NEresp_new[[j]]
  
  nrow(temp.pred)
  length(temp.resp)
  
  n <- nrow(temp.pred)
  set.seed(156)
  folds <- sample(rep(1:k, length.out = n))
  # Store errors
  mse <- numeric(k)
  for (i in 1:k) {
    # Split into training and test sets
    test.id <- which(folds == i)
    train.preds <- temp.pred[-test.id, ]
    train.resp <- temp.resp[-test.id]
    test.preds <- temp.pred[test.id, ]
    test.resp <- temp.resp[test.id]
    
    # Fit model
    train.data.fit <- as.data.frame(cbind(train.resp, train.preds))
    names(train.data.fit)[1] <- "co"
    
    temp.model <- lm(formula(NErefits[[j]]), train.data.fit)
    
    # Predict
    preds <- predict(temp.model, newdata = test.preds)
    
    # Compute Mean Squared Error
    mse[i] <- mean((test.resp - preds)^2)
  }
  NE.kcv[[j]] <- mean(mse)
}


#SE aus
SE.kcv <- NULL
for(j in 1:3){ 
  temp.pred <- SEpreds_new[[j]]
  temp.resp <- SEresp_new[[j]]
  
  
  nrow(temp.pred)
  length(temp.resp)
  
  n <- nrow(temp.pred)
  set.seed(157)
  folds <- sample(rep(1:k, length.out = n))
  # Store errors
  mse <- numeric(k)
  for (i in 1:k) {
    # Split into training and test sets
    test.id <- which(folds == i)
    train.preds <- temp.pred[-test.id, ]
    train.resp <- temp.resp[-test.id]
    test.preds <- temp.pred[test.id, ]
    test.resp <- temp.resp[test.id]
    
    # Fit model
    train.data.fit <- as.data.frame(cbind(train.resp, train.preds))
    names(train.data.fit)[1] <- "co"
    
    temp.model <- lm(formula(SErefits[[j]]), train.data.fit)
    
    # Predict
    preds <- predict(temp.model, newdata = test.preds)
    
    # Compute Mean Squared Error
    mse[i] <- mean((test.resp - preds)^2)
  }
  SE.kcv[[j]] <- mean(mse)
}

NE.kcv
SE.kcv


#eBIC k-fold
NE.kcv.ebic <- NULL
for (l in 1:length(gamma.seq)) {
  NE.kcv.temp <- NULL
  for(j in 1:3){ 
    temp.pred <- NEpreds_new[[j]]
    temp.resp <- NEresp_new[[j]]
    
    nrow(temp.pred)
    length(temp.resp)
    
    n <- nrow(temp.pred)
    set.seed(156)
    folds <- sample(rep(1:k, length.out = n))
    # Store errors
    mse <- numeric(k)
    for (i in 1:k) {
      # Split into training and test sets
      test.id <- which(folds == i)
      train.preds <- temp.pred[-test.id, ]
      train.resp <- temp.resp[-test.id]
      test.preds <- temp.pred[test.id, ]
      test.resp <- temp.resp[test.id]
      
      # Fit model
      train.data.fit <- as.data.frame(cbind(train.resp, train.preds))
      names(train.data.fit)[1] <- "co"
      
      temp.model <- lm(formula(NErefits.eBIC[[j]][[l]]), train.data.fit)
      
      # Predict
      preds <- predict(temp.model, newdata = test.preds)
      
      # Compute Mean Squared Error
      mse[i] <- mean((test.resp - preds)^2)
    }
    NE.kcv.temp[[j]] <- mean(mse)
  }
  NE.kcv.ebic[[l]] <- NE.kcv.temp
}

SE.kcv.ebic <- NULL
for (l in 1:length(gamma.seq)) {
  SE.kcv.temp <- NULL
  for(j in 1:3){ 
    temp.pred <- SEpreds_new[[j]]
    temp.resp <- SEresp_new[[j]]
    
    
    nrow(temp.pred)
    length(temp.resp)
    
    n <- nrow(temp.pred)
    set.seed(157)
    folds <- sample(rep(1:k, length.out = n))
    # Store errors
    mse <- numeric(k)
    for (i in 1:k) {
      # Split into training and test sets
      test.id <- which(folds == i)
      train.preds <- temp.pred[-test.id, ]
      train.resp <- temp.resp[-test.id]
      test.preds <- temp.pred[test.id, ]
      test.resp <- temp.resp[test.id]
      
      # Fit model
      train.data.fit <- as.data.frame(cbind(train.resp, train.preds))
      names(train.data.fit)[1] <- "co"
      
      temp.model <- lm(formula(SErefits.eBIC[[j]][[l]]), train.data.fit)
      
      # Predict
      preds <- predict(temp.model, newdata = test.preds)
      
      # Compute Mean Squared Error
      mse[i] <- mean((test.resp - preds)^2)
    }
    SE.kcv.temp[[j]] <- mean(mse)
  }
  
  SE.kcv.ebic[[l]] <- SE.kcv.temp
  
}

SE.kcv.ebic[[1]]
SE.kcv

#output/save
setwd("~/CO_AUS/Aus_CO-main/Interactions_New")
save(NE.kcv, NE.kcv.ebic, 
     SE.kcv, SE.kcv.ebic,
     file = "validation_kfold.rda")

