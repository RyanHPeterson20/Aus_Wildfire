##for performing cross-validation work with (fused) lasso

#How do we want to arrange all of this??

##k-fold
##loo-cv
##predictions with above information

#libraries
suppressMessages(library(genlasso)) #used for fused lasso 
suppressMessages(library(MASS)) #for matrix shenanigans in functions.

# Data Set-Up
setwd("~/CO_AUS/Aus_CO-main")

load( "ne_data.rda")
load( "se_data.rda")
load( "bounded_data.rda")
load( "data_matrix.rda")
load( "lag_list.rda")

source("group_functions.R")
source("lasso_valid_functions.R") #include predict.fusedlasso
source("genlasso_func.R") #includes cv.fusedlasso (DO NOT USE in final work, correct and move)

#set up season years/weeks
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

#group response matrices
#NEAus
NEAus_1 <- NEbase_matrix[ ,1:3]
NEAus_2 <- NEbase_matrix[ ,4:8]
NEAus_3 <- NEbase_matrix[ ,9:12]
NEAus_4 <- NEbase_matrix[ ,13:17]
NEAus_5 <- NEbase_matrix[ ,18:21]
NEAus_6 <- NEbase_matrix[ ,22:32]

#SEAus
SEAus_1 <- SEbase_matrix[ ,1:3]
SEAus_2 <- SEbase_matrix[ ,4:7]
SEAus_3 <- SEbase_matrix[ ,8:16]
SEAus_4 <- SEbase_matrix[ ,17:20]
SEAus_5 <- SEbase_matrix[ ,21:32]

NEAus_mat <- list(NEAus_1, NEAus_2, NEAus_3, NEAus_4, NEAus_5, NEAus_6)
SEAus_mat <- list(SEAus_1, SEAus_2, SEAus_3, SEAus_4, SEAus_5)

rm(NEAus_1, NEAus_2, NEAus_3, NEAus_4, NEAus_5, NEAus_6, SEAus_1, SEAus_2, SEAus_3, SEAus_4, SEAus_5)

# distance matrices
D <- getD1d(52)
empty_D <- matrix(rep(0, (52*51)) , ncol = 52, nrow = 51)
D1 <- cbind(D, empty_D, empty_D, empty_D)
D2 <- cbind(empty_D, D, empty_D, empty_D)
D3 <- cbind(empty_D, empty_D, D, empty_D)
D4 <- cbind(empty_D, empty_D, empty_D, D)

D_new <- rbind(D1, D2, D3, D4)

#Update for OLR (D5)
D1 <- cbind(D, empty_D, empty_D, empty_D, empty_D)
D2 <- cbind(empty_D, D, empty_D, empty_D, empty_D)
D3 <- cbind(empty_D, empty_D, D, empty_D, empty_D)
D4 <- cbind(empty_D, empty_D, empty_D, D, empty_D)
D5 <- cbind(empty_D, empty_D, empty_D, empty_D, D)

D_olr <- rbind(D1, D2, D3, D4, D5)

rm(D, empty_D, D1, D2, D3, D4, D5)

# grouping functions 

#full model
NE_preds <- NElag_grouping(NE_laglist = NE_laglist_std, j = -c(19))
NE_resp <- NEresp_grouping(NEAus_mat = NEAus_mat, j = -c(19))

SE_preds <- SElag_grouping(SE_laglist = SE_laglist_std, j = -c(19))
SE_resp <- SEresp_grouping(SEAus_mat = SEAus_mat, j = -c(19))

NE_gamma <- 0.85
NEfuse_grouplist <- list()
NEfuse_cv <- list()
NE_lambdamin <- c()

n <- length(NE_resp)
for (i in 1:n) {
  NEresp_temp <- NE_resp[[i]]
  NEpred_temp <- as.matrix(NE_preds[[i]][ ,1:260])
  
  NEgroup_temp <- fusedlasso(y = NEresp_temp, X = NEpred_temp, D_olr, gamma = NE_gamma)
  
  NEgroup_cv <- cv.fusedlasso(NEgroup_temp, k =5, D_olr) 
    
  NE_lambdamin <- c(NE_lambdamin, NEgroup_cv$lambda.min)  
  NEfuse_grouplist[[paste0("Group_", i)]] <- NEgroup_temp
  NEfuse_cv[[paste0("Group_", i)]] <- NEgroup_cv
}


SE_gamma <- 0.85
SEfuse_grouplist <- list()
SEfuse_cv <- list()
SE_lambdamin <- c()

n <- length(SE_resp)
for (i in 1:n) {
  SEresp_temp <- SE_resp[[i]]
  SEpred_temp <- as.matrix(SE_preds[[i]][ ,1:260])
  
  SEgroup_temp <- fusedlasso(y = SEresp_temp, X = SEpred_temp, D_olr, gamma = SE_gamma)
  
  SEgroup_cv <- cv.fusedlasso(SEgroup_temp, k =5, D_olr) 
  
  SE_lambdamin <- c(SE_lambdamin, SEgroup_cv$lambda.min) 
  SEfuse_grouplist[[paste0("Group_", i)]] <- SEgroup_temp
  SEfuse_cv[[paste0("Group_", i)]] <- SEgroup_cv
}



# compare lambda's

#NE Aus
plot(NEfuse_cv[[1]])
plot(NEfuse_cv[[2]])
plot(NEfuse_cv[[3]])
plot(NEfuse_cv[[4]])
plot(NEfuse_cv[[5]])
plot(NEfuse_cv[[6]])

#SE Aus
plot(SEfuse_cv[[1]])
plot(SEfuse_cv[[2]]) #use one stand error
plot(SEfuse_cv[[3]])
plot(SEfuse_cv[[4]])
plot(SEfuse_cv[[5]])



NE_newlambda <- NE_lambdamin
SE_newlambda <- SE_lambdamin

new_min <- which.min(SEfuse_cv[[1]]$err[1:800])
alt_SEgroup1_lambda <- SEfuse_cv$Group_1$lambda[new_min] #alternative lambda
SE_newlambda[1] <- SEfuse_cv$Group_1$lambda[new_min]

SE_newlambda[2] <- SEfuse_cv[[2]]$lambda.1se

NE_oldlambda <- c()
for(k in 1:length(NEfuse_grouplist)){
  n <- length(NE_resp[[j]])
  p <- ncol(NE_preds[[j]])
  lambda_temp <- sqrt(n * log(p))
  NE_oldlambda <- c(NE_oldlambda, lambda_temp)
}

SE_oldlambda <- c()
for(k in 1:length(SEfuse_grouplist)){
  n <- length(SE_resp[[j]])
  p <- ncol(SE_preds[[j]])
  lambda_temp <- sqrt(n * log(p))
  SE_oldlambda <- c(SE_oldlambda, lambda_temp)
}

#TODO: add legend if we are going to use this.
plot(1:6, NE_newlambda, type = "b", pch = 16, 
     xlab = "Groups", ylab = "lambda", main = "NE Aus Lambda values" )
abline(h = NE_oldlambda[1], lty =2, col = "red2")

plot(1:5, SE_newlambda, type = "b", pch = 16, ylim = range(SE_newlambda, SE_oldlambda),
     xlab = "Groups", ylab = "lambda", main = "SE Aus Lambda values" )
abline(h = SE_oldlambda[1], lty =2, col = "red2")


#TODO: complete loocv for variable selection

NEfuse_all <- list()
SEfuse_all <- list()

for(k in 1:18){
  NE_preds <- NElag_grouping(NE_laglist = NE_laglist_std, j = -c(19,k))
  NE_resp <- NEresp_grouping(NEAus_mat = NEAus_mat, j = -c(19,k))
  
  SE_preds <- SElag_grouping(SE_laglist = SE_laglist_std, j = -c(19,k))
  SE_resp <- SEresp_grouping(SEAus_mat = SEAus_mat, j = -c(19,k))
  
  NE_gamma <- 0.85
  NEfuse_loo <- list()
  
  n <- length(NE_resp)
  for (i in 1:n) {
    NEresp_temp <- NE_resp[[i]]
    NEpred_temp <- as.matrix(NE_preds[[i]][ ,1:208])
    
    NEgroup_temp <- fusedlasso(y = NEresp_temp, X = NEpred_temp, D_new, gamma = NE_gamma)
    
    NEfuse_loo[[paste0("Group_", i)]] <- NEgroup_temp
  }
  
  NEfuse_all[[paste0("TestYear", season_years[k])]] <- NEfuse_loo
  
  SE_gamma <- 0.85
  SEfuse_loo <- list()
  
  n <- length(SE_resp)
  for (i in 1:n) {
    SEresp_temp <- SE_resp[[i]]
    SEpred_temp <- as.matrix(SE_preds[[i]][ ,1:208])
    
    SEgroup_temp <- fusedlasso(y = SEresp_temp, X = SEpred_temp, D_new, gamma = SE_gamma)
    
    SEfuse_loo[[paste0("Group_", i)]] <- SEgroup_temp
  }
  
  SEfuse_all[[paste0("TestYear", season_years[k])]] <- SEfuse_loo
  
}

NEfuse_all_new <- NEfuse_all
NEfuse_all_new <- SEfuse_all

save(NEfuse_all_new, SEfuse_all_new, NE_newlambda, SE_newlambda,
     NEfuse_grouplist, SEfuse_grouplist, file = "lasso_pred.rda")

save(NEfuse_all, SEfuse_all, file = "lasso_loo.rda") #save 

#compare coefficients

#baseline (full model coefs)

NE_coef <- list()
NE_ratio <- list()
for (j in 1:6) {
  NE_beta <- coef(NEfuse_grouplist[[j]], lambda = NE_newlambda[j])$beta
  
  NEgroup_coef <- matrix(0, ncol = length(NE_beta))
  NEcoef_count <- matrix(0, ncol = length(NE_beta)) 
  colnames(NEgroup_coef) <- colnames(NE_preds[[j]][1:208])
  colnames(NEcoef_count) <- colnames(NEgroup_coef)
  
  for (k in 1:18) {
    temp_fit <- NEfuse_all[[k]][[j]]
    temp_beta <- coef(temp_fit, lambda = NE_newlambda[j])$beta
    NEgroup_coef <- rbind(NEgroup_coef, t(temp_beta))
    
    coef_index <- which(temp_beta != 0)
    coeff <- numeric(208)
    coeff[coef_index] <- 1
    NEcoef_count <- rbind(NEcoef_count, coeff)
  }
  
  NEgroup_coef <- NEgroup_coef[-1, ]
  
  NE_coef[[paste0("Group_", j)]] <- NEgroup_coef
  
  coef_ratio <- colSums(NEcoef_count)/18
  NE_ratio[[paste0("Group_", j)]] <- coef_ratio[which(NE_beta != 0)]

  
}

#SE Aus
SE_coef <- list()
SE_ratio <- list()
for (j in 1:5) {
  SE_beta <- coef(SEfuse_grouplist[[j]], lambda = SE_newlambda[j])$beta
  
  SEgroup_coef <- matrix(0, ncol = length(SE_beta))
  SEcoef_count <- matrix(0, ncol = length(SE_beta)) 
  colnames(SEgroup_coef) <- colnames(SE_preds[[j]][1:208])
  colnames(SEcoef_count) <- colnames(SEgroup_coef)
  
  for (k in 1:18) {
    temp_fit <- SEfuse_all[[k]][[j]]
    temp_beta <- coef(temp_fit, lambda = SE_newlambda[j])$beta
    SEgroup_coef <- rbind(SEgroup_coef, t(temp_beta))
    
    coef_index <- which(temp_beta != 0)
    coeff <- numeric(208)
    coeff[coef_index] <- 1
    SEcoef_count <- rbind(SEcoef_count, coeff)
  }
  
  SEgroup_coef <- SEgroup_coef[-1, ]
  
  SE_coef[[paste0("Group_", j)]] <- SEgroup_coef
  
  coef_ratio <- colSums(SEcoef_count)/18
  SE_ratio[[paste0("Group_", j)]]   <- coef_ratio[which(SE_beta != 0)]

}

#TODO: clean up the below plots
boxplot(NE_ratio, ylim = c(0,1), main = "(NE Aus) LOO Coef Ratio")
boxplot(SE_ratio, ylim = c(0,1), main = "(SE Aus) LOO Coef Ratio")

#probably not going to use these
#TODO: add all together for a NE and SE hist.
hist(NE_ratio$Group_1, xlim = c(0,1))  
hist(NE_ratio$Group_2, xlim = c(0,1))  
hist(NE_ratio$Group_3, xlim = c(0,1))
#....
hist(SE_ratio$Group_4, xlim = c(0,1))
hist(SE_ratio$Group_5, xlim = c(0,1))


#get predictions for each test year (loo) with prediction CI (and combined R^2)


#TODO: move predictions over to a different file
#need to save full fused, loo fused lasso, cv lambdas
#import methods of groupings

#include by week predictions instead of by group

#begin with 2019 (work backwards)

#and generalize for out of sample preds 

NE_testpreds <- NElag_grouping(NE_laglist = NE_laglist_std, j = -c(1:18))
NE_testresp <- NEresp_grouping(NEAus_mat = NEAus_mat, j = -c(1:18))

NE_2019preds <- as.data.frame(matrix(NA, ncol = 4))
colnames(NE_2019preds) <- c("y"," y.hat", "PI.upper", "PI.lower")

for (j in 1:6) {
  test_object <- NEfuse_grouplist[[j]]
  test_lambda <- NE_newlambda[j]
  test_resp <- NE_testresp[[j]]
  test_preds <- NE_testpreds[[j]][,1:260]
  
  test_out <- predict.fusedlasso(test_object, test_lambda, 
                                 y_new = test_resp, X_new = test_preds)
  NE_2019preds <- rbind(NE_2019preds, test_out)
}

NE_2019preds <- NE_2019preds[-1, ]

x_vals <- 1:32
NEy_range <- range(NE_2019preds)

plot(x_vals, NE_2019preds$y, type = "l", lwd = 2, ylim = NEy_range)
lines(x_vals, NE_2019preds$` y.hat`, lwd = 2, col = "magenta3")
lines(x_vals, NE_2019preds$PI.upper, lwd = 1.5, lty = 2,
      col = "red2")
lines(x_vals, NE_2019preds$PI.lower, lwd = 1.5, lty = 2,
      col = "blue2")

#SE Aus 2019/2020 preds
SE_testpreds <- SElag_grouping(SE_laglist = SE_laglist_std, j = -c(1:18))
SE_testresp <- SEresp_grouping(SEAus_mat = SEAus_mat, j = -c(1:18))

SE_2019preds <- as.data.frame(matrix(NA, ncol = 4))
colnames(SE_2019preds) <- c("y"," y.hat", "PI.upper", "PI.lower")

for (j in 1:5) {
  test_object <- SEfuse_grouplist[[j]]
  test_lambda <- SE_newlambda[j]
  test_resp <- SE_testresp[[j]]
  test_preds <- SE_testpreds[[j]][,1:260]
  
  test_out <- predict.fusedlasso(test_object, test_lambda, 
                                 y_new = test_resp, X_new = test_preds)
  SE_2019preds <- rbind(SE_2019preds, test_out)
}

SE_2019preds <- SE_2019preds[-1, ]

x_vals <- 1:32
SEy_range <- range(SE_2019preds)

plot(x_vals, SE_2019preds$y, type = "l", lwd = 2, ylim = SEy_range)
lines(x_vals, SE_2019preds$` y.hat`, lwd = 2, col = "magenta3")
lines(x_vals, SE_2019preds$PI.upper, lwd = 1.5, lty = 2,
      col = "red2")
lines(x_vals, SE_2019preds$PI.lower, lwd = 1.5, lty = 2,
      col = "blue2")




#test sections
j <- 2 #group
test_object <- NEfuse_grouplist[[j]]
test_lambda <- NE_newlambda[j]
test_resp <- NE_testresp[[j]]
test_preds <- NE_testpreds[[j]][,1:208]

test_out <- predict.fusedlasso(test_object, test_lambda, 
                               y_new = test_resp, X_new = test_preds)

x_vals <- 1:nrow(test_out)
y_range <- range(test_out)

plot(x_vals, test_out$y, type = "l", lwd = 2, ylim = y_range)
lines(x_vals, test_out$y.hat, lwd = 2, col = "magenta3")
lines(x_vals, test_out$PI.upper, lwd = 1.5, lty = 2,
      col = "red2")
lines(x_vals, test_out$PI.lower, lwd = 1.5, lty = 2,
      col = "blue2")



#testing preds
years <- 1:18
k <- years[-2]

NE_testpreds <- NElag_grouping(NE_laglist = NE_laglist_std, j = -c(19, k))
NE_testresp <- NEresp_grouping(NEAus_mat = NEAus_mat, j = -c(19, k))

j <- 1 #group
temp_coef <- NE_coef[[j]][2, ]
temp_cov <- as.matrix(NE_testpreds[[j]][,1:208])
temp_cov %*% temp_coef
NE_testresp[[j]]

#TODO: modifiy cv.fusedlasso test fold predictions and R^2


#intercept vector (set-up) Temporary
#TODO: testing a manual intercept term
## int_vector <- rep(1, length.out = length(NE_preds[[1]]$nino_lag1))
## NE_preds_int <- data.frame(int = int_vector, NE_preds[[1]])

## D_int <- cbind(D_new, rep(0, length.out = 204)) 
## D_int <- rbind( c(1, rep(0, length.out = 208)), D_int)


#TODO: compare the different version of lambda
#e.g. the version the is empirically defined and mse min
#below is emp lambda coef work: 
j <- 1
#for(k in 1:length(NEfuse_grouplist)){
n <- length(NE_resp[[j]])
p <- ncol(NE_preds[[j]])
lambda1 <- sqrt(n * log(p))

NEfuse_grouplist[[j]]$lambda

NEgroup_coef <- coef(NEfuse_grouplist[[j]], lambda = lambda1)
NEfull_coef <- which(NEgroup_coef$beta != 0)



#TODO: double check with various version of lambda 
j <- 1
SEgroup_coef <- coef(SEfuse_grouplist[[j]], lambda = SE_newlambda[j])$beta
SEgroup1_nino <- SEgroup_coef[1:52]
SEgroup1_dmi <- SEgroup_coef[53:104]
SEgroup1_tsa <- SEgroup_coef[105:156]
SEgroup1_aao <- SEgroup_coef[157:208]



#TEST plots for new lambdas
#SE group 1
#png(filename = "SE_group1.png", width = 2400, height = 1800, res = 300)
par(mar = c(5, 4, 4, 9) + 0.2, xpd = TRUE)  # Increase the right margin (fourth value) to make space for the legend
plot(1:52, SEgroup1_nino, type = "l", xlab = "Lag", ylab = "Coefficients",
     main = "SE Group 1", ylim = range(SEgroup_coef), 
     col = "darkmagenta")
#lines(1:52, SE_coefs_base$Group_1$nino_base, col = "magenta", lty = 2, lwd = 1.5)
lines(1:52, SEgroup1_dmi, col = "blue")
#lines(1:52, SE_coefs_base$Group_1$dmi_base, col = "cyan", lty = 2, lwd = 1.5)
lines(1:52, SEgroup1_tsa, col = "red")
#lines(1:52, SE_coefs_base$Group_1$tsa_base, col = "coral", lty = 2, lwd = 1.5)
lines(1:52, SEgroup1_aao, col = "darkgreen")
#lines(1:52, SE_coefs_base$Group_1$aao_base, col = "green", lty = 2, lwd = 1.5)
segments(x0 = -1, y0 = 0, x1 = 54, y1 = 0, lty = 2)
legend("topright", inset = c(-0.3, 0),
       legend = legend1_names,
       col = legend1_col, 
       lty = legend1_lty, lwd = 2, bty = "o") #, bg = "white")
dev.off()




#begin with 2018 (work backwards)
k <- c(1:17,19)

NE_testpreds <- NElag_grouping(NE_laglist = NE_laglist_std, j = -k)
NE_testresp <- NEresp_grouping(NEAus_mat = NEAus_mat, j = -k)

SE_testpreds <- SElag_grouping(SE_laglist = SE_laglist_std, j = -k)
SE_testresp <- SEresp_grouping(SEAus_mat = SEAus_mat, j = -k)

NE_2018preds <- as.data.frame(matrix(NA, ncol = 4))
colnames(NE_2018preds) <- c("y"," y.hat", "PI.upper", "PI.lower")

for (j in 1:6) {
  test_object <- NEfuse_grouplist[[j]]
  test_lambda <- NE_newlambda[j]
  test_resp <- NE_testresp[[j]]
  test_preds <- NE_testpreds[[j]][,1:208]
  
  test_out <- predict.fusedlasso(test_object, test_lambda, 
                                 y_new = test_resp, X_new = test_preds)
  
  NE_2018preds <- rbind(NE_2018preds, test_out)
}

NE_2018preds <- NE_2018preds[-1, ]

x_vals <- 1:32
NEy_range <- range(NE_2018preds)

plot(x_vals, NE_2018preds$y, type = "l", lwd = 2, ylim = NEy_range)
lines(x_vals, NE_2018preds$` y.hat`, lwd = 2, col = "magenta3")
lines(x_vals, NE_2018preds$PI.upper, lwd = 1.5, lty = 2,
      col = "red2")
lines(x_vals, NE_2018preds$PI.lower, lwd = 1.5, lty = 2,
      col = "blue2")

#SE Aus 2019/2020 preds


SE_2018preds <- as.data.frame(matrix(NA, ncol = 4))
colnames(SE_2018preds) <- c("y"," y.hat", "PI.upper", "PI.lower")

for (j in 1:5) {
  test_object <- SEfuse_grouplist[[j]]
  test_lambda <- SE_newlambda[j]
  test_resp <- SE_testresp[[j]]
  test_preds <- SE_testpreds[[j]][,1:208]
  
  test_out <- predict.fusedlasso(test_object, test_lambda, 
                                 y_new = test_resp, X_new = test_preds)
  SE_2018preds <- rbind(SE_2018preds, test_out)
}

SE_2018preds <- SE_2018preds[-1, ]

x_vals <- 1:32
SEy_range <- range(SE_2018preds)

plot(x_vals, SE_2018preds$y, type = "l", lwd = 2, ylim = SEy_range)
lines(x_vals, SE_2018preds$` y.hat`, lwd = 2, col = "magenta3")
lines(x_vals, SE_2018preds$PI.upper, lwd = 1.5, lty = 2,
      col = "red2")
lines(x_vals, SE_2018preds$PI.lower, lwd = 1.5, lty = 2,
      col = "blue2")




##--------------------quantile section--------------------##


#TODO: work include quantile data

setwd("~/CO_AUS/Aus_CO-main")
load( "data_quantile.rda") #quantile data

#first include indicator quantiles as series for each climate mode.
# distance matrices 

#TODO: update distance matrices for indicator quantiles $I[x>c]$
##using D for nino, nino_q, dmi, dmi_q

D <- getD1d(52)
empty_D <- matrix(rep(0, (52*51)) , ncol = 52, nrow = 51)
D1 <- cbind(D, empty_D, empty_D, empty_D) #nino
D2 <- cbind(empty_D, D, empty_D, empty_D) #nino_q
D3 <- cbind(empty_D, empty_D, D, empty_D) #dmi
D4 <- cbind(empty_D, empty_D, empty_D, D) #dmi_q

D_new <- rbind(D1, D2, D3, D4)

#TODO: change for other variations later (including simple olr)
#Update for OLR (D5) 
D1 <- cbind(D, empty_D, empty_D, empty_D, empty_D)
D2 <- cbind(empty_D, D, empty_D, empty_D, empty_D)
D3 <- cbind(empty_D, empty_D, D, empty_D, empty_D)
D4 <- cbind(empty_D, empty_D, empty_D, D, empty_D)
D5 <- cbind(empty_D, empty_D, empty_D, empty_D, D)

D_olr <- rbind(D1, D2, D3, D4, D5)

rm(D, empty_D, D1, D2, D3, D4, D5)

# grouping functions 

#full model
NE_preds <- NElag_grouping(NE_laglist = NE_laglist_std, j = -c(19))
NE_resp <- NEresp_grouping(NEAus_mat = NEAus_mat, j = -c(19))

SE_preds <- SElag_grouping(SE_laglist = SE_laglist_std, j = -c(19))
SE_resp <- SEresp_grouping(SEAus_mat = SEAus_mat, j = -c(19))

#quantile 90 indicator
NE_preds_q90 <- NElag_grouping(NE_laglist = NE_laglist_q90, j = -c(19))
SE_preds_q90 <- SElag_grouping(SE_laglist = SE_laglist_q90, j = -c(19))

#quantile 75 indicator
NE_preds_q75 <- NElag_grouping(NE_laglist = NE_laglist_q75, j = -c(19))
SE_preds_q75 <- SElag_grouping(SE_laglist = SE_laglist_q75, j = -c(19))

#cv checks

NE_gamma <- 0.85
NEfuse_grouplist <- list()
NEfuse_cv <- list()
NE_lambdamin <- c()

n <- length(NE_resp)
for (i in 1:n) {
  NEresp_temp <- NE_resp[[i]]
  NEpred_temp <- as.matrix(cbind(NE_preds[[i]][ ,1:52], NE_preds_q75[[i]][ ,1:52],
                                 NE_preds[[i]][ ,53:104], NE_preds_q75[[i]][ ,53:104]))
  
  NEgroup_temp <- fusedlasso(y = NEresp_temp, X = NEpred_temp, D_new, gamma = NE_gamma)
  
  NEgroup_cv <- cv.fusedlasso(NEgroup_temp, k =5, D_new) 
  
  NE_lambdamin <- c(NE_lambdamin, NEgroup_cv$lambda.min)  
  NEfuse_grouplist[[paste0("Group_", i)]] <- NEgroup_temp
  NEfuse_cv[[paste0("Group_", i)]] <- NEgroup_cv
}

SE_gamma <- 0.85
SEfuse_grouplist <- list()
SEfuse_cv <- list()
SE_lambdamin <- c()

n <- length(SE_resp)
for (i in 1:n) {
  SEresp_temp <- SE_resp[[i]]
  SEpred_temp <- as.matrix(cbind(SE_preds[[i]][ ,1:52], SE_preds_q75[[i]][ ,1:52],
                                 SE_preds[[i]][ ,53:104], SE_preds_q75[[i]][ ,53:104]))
  
  SEgroup_temp <- fusedlasso(y = SEresp_temp, X = SEpred_temp, D_new, gamma = SE_gamma)
  
  SEgroup_cv <- cv.fusedlasso(SEgroup_temp, k =5, D_new) 
  
  SE_lambdamin <- c(SE_lambdamin, SEgroup_cv$lambda.min) 
  SEfuse_grouplist[[paste0("Group_", i)]] <- SEgroup_temp
  SEfuse_cv[[paste0("Group_", i)]] <- SEgroup_cv
}

#lambda checks

#NE Aus
plot(NEfuse_cv[[1]])
plot(NEfuse_cv[[2]])
plot(NEfuse_cv[[3]])
plot(NEfuse_cv[[4]])
plot(NEfuse_cv[[5]])
plot(NEfuse_cv[[6]])

#SE Aus
plot(SEfuse_cv[[1]])
plot(SEfuse_cv[[2]]) 
plot(SEfuse_cv[[3]])
plot(SEfuse_cv[[4]])
plot(SEfuse_cv[[5]])

NE_newlambda <- NE_lambdamin
SE_newlambda <- SE_lambdamin

#coef viz (from AGU_viz)

#NE Visualizations
NE_coefs <- list()
#full coef list
NEfuse_coefs <- list()
NE_range <- c()

for(k in 1:length(NEfuse_grouplist)){
  
  lambda1 <- NE_newlambda[k]
  
  NEgroup_coef <- coef(NEfuse_grouplist[[k]], lambda = lambda1)
  
  NEfuse <- NEgroup_coef$beta
  colnames(NEfuse) <- c("NEgroup")
  NEfuse_coefs[[paste0("Group_", k)]] <- as.data.frame(NEfuse)
  
  NE_range <- c(NE_range, range(NEgroup_coef$beta))
  
  #extract each index
  nino_coef <- NEgroup_coef$beta[1:52,]
  ninoQ_coef <-    NEgroup_coef$beta[53:104,]
  dmi_coef <- NEgroup_coef$beta[105:156,]
  dmiQ_coef <- NEgroup_coef$beta[157:208,]
  
  NE_coefs[[paste0("Group_", k)]] <- data.frame(nino_coef, ninoQ_coef, dmi_coef, dmiQ_coef)
}

#plot ne coefs

#group1
plot(1:52, NE_coefs$Group_1$nino_coef, type = "l", xlab = "Lag", ylab = "Coefficients",
     main = "NE Group 1", ylim = range(NE_range))
lines(1:52, NE_coefs$Group_1$ninoQ_coef, col = "blue")
lines(1:52, NE_coefs$Group_1$dmi_coef, col = "red")
lines(1:52, NE_coefs$Group_1$dmiQ_coef, col = "darkgreen")

#group 2
plot(1:52, NE_coefs$Group_2$nino_coef, type = "l", xlab = "Lag", ylab = "Coefficients",
     main = "NE Group 2", ylim = range(NE_range))
lines(1:52, NE_coefs$Group_2$ninoQ_coef, col = "blue")
lines(1:52, NE_coefs$Group_2$dmi_coef, col = "red")
lines(1:52, NE_coefs$Group_2$dmiQ_coef, col = "darkgreen")

#group 3
plot(1:52, NE_coefs$Group_3$nino_coef, type = "l", xlab = "Lag", ylab = "Coefficients",
     main = "NE Group 3", ylim = range(NE_range))
lines(1:52, NE_coefs$Group_3$ninoQ_coef, col = "blue")
lines(1:52, NE_coefs$Group_3$dmi_coef, col = "red")
lines(1:52, NE_coefs$Group_3$dmiQ_coef, col = "darkgreen")

#group 4
plot(1:52, NE_coefs$Group_4$nino_coef, type = "l", xlab = "Lag", ylab = "Coefficients",
     main = "NE Group 4", ylim = range(NE_range))
lines(1:52, NE_coefs$Group_4$ninoQ_coef, col = "blue")
lines(1:52, NE_coefs$Group_4$dmi_coef, col = "red")
lines(1:52, NE_coefs$Group_4$dmiQ_coef, col = "darkgreen")

#group 5
plot(1:52, NE_coefs$Group_5$nino_coef, type = "l", xlab = "Lag", ylab = "Coefficients",
     main = "NE Group 5", ylim = range(NE_range))
lines(1:52, NE_coefs$Group_5$ninoQ_coef, col = "blue")
lines(1:52, NE_coefs$Group_5$dmi_coef, col = "red")
lines(1:52, NE_coefs$Group_5$dmiQ_coef, col = "darkgreen")

#group 6
plot(1:52, NE_coefs$Group_6$nino_coef, type = "l", xlab = "Lag", ylab = "Coefficients",
     main = "NE Group 6", ylim = range(NE_range))
lines(1:52, NE_coefs$Group_6$ninoQ_coef, col = "blue")
lines(1:52, NE_coefs$Group_6$dmi_coef, col = "red")
lines(1:52, NE_coefs$Group_6$dmiQ_coef, col = "darkgreen")


#se coefs
#SE Visualizations
SE_coefs <- list()
#full coef list
SEfuse_coefs <- list()

SE_range <- c()

for(k in 1:length(SEfuse_grouplist)){
  
  lambda1 <- SE_newlambda[k]
  
  SEgroup_coef <- coef(SEfuse_grouplist[[k]], lambda = lambda1)
  
  SEfuse <- SEgroup_coef$beta
  colnames(SEfuse) <- c("SEgroup")
  SEfuse_coefs[[paste0("Group_", k)]] <- as.data.frame(SEfuse)
  
  SE_range <- c(SE_range, range(SEgroup_coef$beta))
  

  #extract each index
  nino_coef <- SEgroup_coef$beta[1:52,]
  ninoQ_coef <-    SEgroup_coef$beta[53:104,]
  dmi_coef <- SEgroup_coef$beta[105:156,]
  dmiQ_coef <- SEgroup_coef$beta[157:208,]
  
  SE_coefs[[paste0("Group_", k)]] <- data.frame(nino_coef, ninoQ_coef, dmi_coef, dmiQ_coef)
  
}



#group1
plot(1:52, SE_coefs$Group_1$nino_coef, type = "l", xlab = "Lag", ylab = "Coefficients",
     main = "SE Group 1", ylim = range(SE_range))
lines(1:52, SE_coefs$Group_1$ninoQ_coef, col = "blue")
lines(1:52, SE_coefs$Group_1$dmi_coef, col = "red")
lines(1:52, SE_coefs$Group_1$dmiQ_coef, col = "darkgreen")


#group 2
plot(1:52, SE_coefs$Group_2$nino_coef, type = "l", xlab = "Lag", ylab = "Coefficients",
     main = "SE Group 2", ylim = range(SE_range))
lines(1:52, SE_coefs$Group_2$ninoQ_coef, col = "blue")
lines(1:52, SE_coefs$Group_2$dmi_coef, col = "red")
lines(1:52, SE_coefs$Group_2$dmiQ_coef, col = "darkgreen")

#group 3
plot(1:52, SE_coefs$Group_3$nino_coef, type = "l", xlab = "Lag", ylab = "Coefficients",
     main = "SE Group 3", ylim = range(SE_range))
lines(1:52, SE_coefs$Group_3$ninoQ_coef, col = "blue")
lines(1:52, SE_coefs$Group_3$dmi_coef, col = "red")
lines(1:52, SE_coefs$Group_3$dmiQ_coef, col = "darkgreen")

#group 4
plot(1:52, SE_coefs$Group_4$nino_coef, type = "l", xlab = "Lag", ylab = "Coefficients",
     main = "SE Group 4", ylim = range(SE_range))
lines(1:52, SE_coefs$Group_4$ninoQ_coef, col = "blue")
lines(1:52, SE_coefs$Group_4$dmi_coef, col = "red")
lines(1:52, SE_coefs$Group_4$dmiQ_coef, col = "darkgreen")

#group 5
plot(1:52, SE_coefs$Group_5$nino_coef, type = "l", xlab = "Lag", ylab = "Coefficients",
     main = "SE Group 5", ylim = range(SE_range))
lines(1:52, SE_coefs$Group_5$ninoQ_coef, col = "blue")
lines(1:52, SE_coefs$Group_5$dmi_coef, col = "red")
lines(1:52, SE_coefs$Group_5$dmiQ_coef, col = "darkgreen")

##second connect an edge of a quantile to its lag.



### --------- modified penalty matrix --------------###

#TODO: rebuild D matrix with different edge connections
#repeat for better connections (start with a single mode (nino) first)
#D_test <- getD1d(52)

D <- getD1d(52)
empty_D <- matrix(rep(0, (52*51)) , ncol = 52, nrow = 51)
D1 <- cbind(D, empty_D, empty_D, empty_D, empty_D)
D2 <- cbind(empty_D, D, empty_D, empty_D, empty_D)
D3 <- cbind(empty_D, empty_D, D, empty_D, empty_D)
D4 <- cbind(empty_D, empty_D, empty_D, D, empty_D)
D5 <- cbind(empty_D, empty_D, empty_D, empty_D, D)

D_olr <- rbind(D1, D2, D3, D4, D5)

rm(D, empty_D, D1, D2, D3, D4, D5)

graph_obj <- getGraph(D_olr)

#TODO: igraph stuff?? (yes)
#testing the inclusion of new nodes and edges
new_graph <- add_vertices(graph_obj, 104)

#create list of alternating values
v1 <- 1:104
v2 <- 261:364
v_c <- c(rbind(v1,v2))

new_graph <- add_edges(new_graph, v_c)

plot(new_graph)

#back to the penalty matrix
D_new <- getDg(new_graph)

#test for nino and indicator nino
#full model
NE_preds <- NElag_grouping(NE_laglist = NE_laglist_std, j = -c(19))
NE_resp <- NEresp_grouping(NEAus_mat = NEAus_mat, j = -c(19))

SE_preds <- SElag_grouping(SE_laglist = SE_laglist_std, j = -c(19))
SE_resp <- SEresp_grouping(SEAus_mat = SEAus_mat, j = -c(19))

#quantile 90 indicator
NE_preds_q90 <- NElag_grouping(NE_laglist = NE_laglist_q90, j = -c(19))
SE_preds_q90 <- SElag_grouping(SE_laglist = SE_laglist_q90, j = -c(19))

#quantile 75 indicator
NE_preds_q75 <- NElag_grouping(NE_laglist = NE_laglist_q75, j = -c(19))
SE_preds_q75 <- SElag_grouping(SE_laglist = SE_laglist_q75, j = -c(19))

#cv checks

NE_gamma <- 0.85
NEfuse_grouplist <- list()
NEfuse_cv <- list()
NE_lambdamin <- c()

n <- length(NE_resp)
for (i in 1:n) {
  NEresp_temp <- NE_resp[[i]]
  NEpred_temp <- as.matrix(cbind(NE_preds[[i]], NE_preds_q75[[i]][ ,1:104]))
  
  NEgroup_temp <- fusedlasso(y = NEresp_temp, X = NEpred_temp, D_new, gamma = NE_gamma)
  
  NEgroup_cv <- cv.fusedlasso(NEgroup_temp, k =5, D_new) 
  
  NE_lambdamin <- c(NE_lambdamin, NEgroup_cv$lambda.min)  
  NEfuse_grouplist[[paste0("Group_", i)]] <- NEgroup_temp
  NEfuse_cv[[paste0("Group_", i)]] <- NEgroup_cv
}


SE_gamma <- 0.85
SEfuse_grouplist <- list()
SEfuse_cv <- list()
SE_lambdamin <- c()

n <- length(SE_resp)
for (i in 1:n) {
  SEresp_temp <- SE_resp[[i]]
  SEpred_temp <- as.matrix(cbind(SE_preds[[i]], SE_preds_q75[[i]][ ,1:104]))
  
  SEgroup_temp <- fusedlasso(y = SEresp_temp, X = SEpred_temp, D_new, gamma = SE_gamma)
  
  SEgroup_cv <- cv.fusedlasso(SEgroup_temp, k =5, D_new) 
  
  SE_lambdamin <- c(SE_lambdamin, SEgroup_cv$lambda.min) 
  SEfuse_grouplist[[paste0("Group_", i)]] <- SEgroup_temp
  SEfuse_cv[[paste0("Group_", i)]] <- SEgroup_cv
}

#lambda checks

#NE Aus
plot(NEfuse_cv[[1]])
plot(NEfuse_cv[[2]])
plot(NEfuse_cv[[3]])
plot(NEfuse_cv[[4]])
plot(NEfuse_cv[[5]])
plot(NEfuse_cv[[6]])

#SE Aus
plot(SEfuse_cv[[1]])
plot(SEfuse_cv[[2]]) 
plot(SEfuse_cv[[3]])
plot(SEfuse_cv[[4]])
plot(SEfuse_cv[[5]])

NE_newlambda <- NE_lambdamin
SE_newlambda <- SE_lambdamin


#NE Visualizations
NE_coefs <- list()
#full coef list
NEfuse_coefs <- list()
NE_range <- c()

for(k in 1:length(NEfuse_grouplist)){
  
  lambda1 <- NE_newlambda[k]
  
  NEgroup_coef <- coef(NEfuse_grouplist[[k]], lambda = lambda1)
  
  NEfuse <- NEgroup_coef$beta
  colnames(NEfuse) <- c("NEgroup")
  NEfuse_coefs[[paste0("Group_", k)]] <- as.data.frame(NEfuse)
  
  NE_range <- c(NE_range, range(NEgroup_coef$beta))
  
  #extract each index
  nino_coef <- NEgroup_coef$beta[1:52,]
  dmi_coef <-    NEgroup_coef$beta[53:104,]
  ninoQ_coef <- NEgroup_coef$beta[105:156,]
  dmiQ_coef <- NEgroup_coef$beta[157:208,]
  
  NE_coefs[[paste0("Group_", k)]] <- data.frame(nino_coef, dmi_coef, ninoQ_coef, dmiQ_coef)
}

#se coefs
#SE Visualizations
SE_coefs <- list()
#full coef list
SEfuse_coefs <- list()

SE_range <- c()

for(k in 1:length(SEfuse_grouplist)){
  
  lambda1 <- SE_newlambda[k]
  
  SEgroup_coef <- coef(SEfuse_grouplist[[k]], lambda = lambda1)
  
  SEfuse <- SEgroup_coef$beta
  colnames(SEfuse) <- c("SEgroup")
  SEfuse_coefs[[paste0("Group_", k)]] <- as.data.frame(SEfuse)
  
  SE_range <- c(SE_range, range(SEgroup_coef$beta))
  
  #extract each index
  nino_coef <- SEgroup_coef$beta[1:52,]
  dmi_coef <-    SEgroup_coef$beta[53:104,]
  ninoQ_coef <- SEgroup_coef$beta[105:156,]
  dmiQ_coef <- SEgroup_coef$beta[157:208,]
  
  SE_coefs[[paste0("Group_", k)]] <- data.frame(nino_coef, dmi_coef, ninoQ_coef, dmiQ_coef)
}

#plot coefs

#ne coefs
#group1
plot(1:52, NE_coefs$Group_1$nino_coef, type = "l", xlab = "Lag", ylab = "Coefficients",
     main = "NE Group 1", ylim = range(NE_range))
abline(h = 0, lty = 2)
lines(1:52, NE_coefs$Group_1$ninoQ_coef, col = "blue")
lines(1:52, NE_coefs$Group_1$dmi_coef, col = "red")
lines(1:52, NE_coefs$Group_1$dmiQ_coef, col = "darkgreen")

#group 2
plot(1:52, NE_coefs$Group_2$nino_coef, type = "l", xlab = "Lag", ylab = "Coefficients",
     main = "NE Group 2", ylim = range(NE_range))
abline(h = 0, lty = 2)
lines(1:52, NE_coefs$Group_2$ninoQ_coef, col = "blue")
lines(1:52, NE_coefs$Group_2$dmi_coef, col = "red")
lines(1:52, NE_coefs$Group_2$dmiQ_coef, col = "darkgreen")

#group 3
plot(1:52, NE_coefs$Group_3$nino_coef, type = "l", xlab = "Lag", ylab = "Coefficients",
     main = "NE Group 3", ylim = range(NE_range))
abline(h = 0, lty = 2)
lines(1:52, NE_coefs$Group_3$ninoQ_coef, col = "blue")
lines(1:52, NE_coefs$Group_3$dmi_coef, col = "red")
lines(1:52, NE_coefs$Group_3$dmiQ_coef, col = "darkgreen")

#group 4
plot(1:52, NE_coefs$Group_4$nino_coef, type = "l", xlab = "Lag", ylab = "Coefficients",
     main = "NE Group 4", ylim = range(NE_range))
abline(h = 0, lty = 2)
lines(1:52, NE_coefs$Group_4$ninoQ_coef, col = "blue")
lines(1:52, NE_coefs$Group_4$dmi_coef, col = "red")
lines(1:52, NE_coefs$Group_4$dmiQ_coef, col = "darkgreen")

#group 5
plot(1:52, NE_coefs$Group_5$nino_coef, type = "l", xlab = "Lag", ylab = "Coefficients",
     main = "NE Group 5", ylim = range(NE_range))
abline(h = 0, lty = 2)
lines(1:52, NE_coefs$Group_5$ninoQ_coef, col = "blue")
lines(1:52, NE_coefs$Group_5$dmi_coef, col = "red")
lines(1:52, NE_coefs$Group_5$dmiQ_coef, col = "darkgreen")

#group 6
plot(1:52, NE_coefs$Group_6$nino_coef, type = "l", xlab = "Lag", ylab = "Coefficients",
     main = "NE Group 6", ylim = range(NE_range))
abline(h = 0, lty = 2)
lines(1:52, NE_coefs$Group_6$ninoQ_coef, col = "blue")
lines(1:52, NE_coefs$Group_6$dmi_coef, col = "red")
lines(1:52, NE_coefs$Group_6$dmiQ_coef, col = "darkgreen")

#se coefs
#group1
plot(1:52, SE_coefs$Group_1$nino_coef, type = "l", xlab = "Lag", ylab = "Coefficients",
     main = "SE Group 1", ylim = range(SE_range))
abline(h = 0, lty = 2)
lines(1:52, SE_coefs$Group_1$ninoQ_coef, col = "blue")
lines(1:52, SE_coefs$Group_1$dmi_coef, col = "red")
lines(1:52, SE_coefs$Group_1$dmiQ_coef, col = "darkgreen")

#group 2
plot(1:52, SE_coefs$Group_2$nino_coef, type = "l", xlab = "Lag", ylab = "Coefficients",
     main = "SE Group 2", ylim = range(SE_range))
abline(h = 0, lty = 2)
lines(1:52, SE_coefs$Group_2$ninoQ_coef, col = "blue")
lines(1:52, SE_coefs$Group_2$dmi_coef, col = "red")
lines(1:52, SE_coefs$Group_2$dmiQ_coef, col = "darkgreen")

#group 3
plot(1:52, SE_coefs$Group_3$nino_coef, type = "l", xlab = "Lag", ylab = "Coefficients",
     main = "SE Group 3", ylim = range(SE_range))
abline(h = 0, lty = 2)
lines(1:52, SE_coefs$Group_3$ninoQ_coef, col = "blue")
lines(1:52, SE_coefs$Group_3$dmi_coef, col = "red")
lines(1:52, SE_coefs$Group_3$dmiQ_coef, col = "darkgreen")

#group 4
plot(1:52, SE_coefs$Group_4$nino_coef, type = "l", xlab = "Lag", ylab = "Coefficients",
     main = "SE Group 4", ylim = range(SE_range))
abline(h = 0, lty = 2)
lines(1:52, SE_coefs$Group_4$ninoQ_coef, col = "blue")
lines(1:52, SE_coefs$Group_4$dmi_coef, col = "red")
lines(1:52, SE_coefs$Group_4$dmiQ_coef, col = "darkgreen")

#group 5
plot(1:52, SE_coefs$Group_5$nino_coef, type = "l", xlab = "Lag", ylab = "Coefficients",
     main = "SE Group 5", ylim = range(SE_range))
abline(h = 0, lty = 2)
lines(1:52, SE_coefs$Group_5$ninoQ_coef, col = "blue")
lines(1:52, SE_coefs$Group_5$dmi_coef, col = "red")
lines(1:52, SE_coefs$Group_5$dmiQ_coef, col = "darkgreen")

