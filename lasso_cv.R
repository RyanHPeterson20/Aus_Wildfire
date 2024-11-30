##for performing cross-validation work with (fused) lasso

#How do we want to arrange all of this??

#libraries
suppressMessages(library(genlasso)) #used for fused lasso 

# Data Set-Up
setwd("~/CO_AUS/Aus_CO-main")

load( "ne_data.rda")
load( "se_data.rda")
load( "bounded_data.rda")
load( "data_matrix.rda")
load( "lag_list.rda")

source("group_functions.R")
source("lasso_valid_functions.R") #include cv.fusedlasso

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

rm(D, empty_D, D1, D2, D3, D4)

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
  NEpred_temp <- as.matrix(NE_preds[[i]][ ,1:208])
  
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
  SEpred_temp <- as.matrix(SE_preds[[i]][ ,1:208])
  
  SEgroup_temp <- fusedlasso(y = SEresp_temp, X = SEpred_temp, D_new, gamma = SE_gamma)
  
  SEgroup_cv <- cv.fusedlasso(SEgroup_temp, k =5, D_new) 
  
  SE_lambdamin <- c(SE_lambdamin, SEgroup_cv$lambda.min) 
  SEfuse_grouplist[[paste0("Group_", i)]] <- SEgroup_temp
  SEfuse_cv[[paste0("Group_", i)]] <- SEgroup_cv
}



# compare lambda's
NE_newlambda <- NE_lambdamin
SE_newlambda <- SE_lambdamin

new_min <- which.min(SEfuse_cv[[1]]$err[1:800])
alt_SEgroup1_lambda <- SEfuse_cv$Group_1$lambda[new_min] #alternative lambda
SE_newlambda[1] <- SEfuse_cv$Group_1$lambda[new_min]

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

save(NEfuse_all, SEfuse_all, file = "lasso_loo.rda") #save 

#compare coefficients

#baseline (full model coefs)
j <- 1 #group number

#for (j in 1:6) {
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
  
  coef_ratio <- colSums(NEcoef_count)/18
  coef_ratio[which(NE_beta != 0)]

#}

#SE Aus
#for (j in 1:5) {
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
  
  coef_ratio <- colSums(SEcoef_count)/18
  coef_ratio[which(SE_beta != 0)]

#}



  
  

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
SEgroup_coef <- coef(SEfuse_grouplist[[1]], lambda = alt_SEgroup1_lambda)

#TEST plots for new lambdas
#SE group 1
png(filename = "SE_group1.png", width = 2400, height = 1800, res = 300)
par(mar = c(5, 4, 4, 9) + 0.2, xpd = TRUE)  # Increase the right margin (fourth value) to make space for the legend
plot(1:52, SE_coefs$Group_1$nino_coef, type = "l", xlab = "Lag", ylab = "Coefficients",
     main = "SE Group 1", ylim = range(SE_range, SEbase_range), 
     col = "darkmagenta")
#lines(1:52, SE_coefs_base$Group_1$nino_base, col = "magenta", lty = 2, lwd = 1.5)
lines(1:52, SE_coefs$Group_1$dmi_coef, col = "blue")
#lines(1:52, SE_coefs_base$Group_1$dmi_base, col = "cyan", lty = 2, lwd = 1.5)
lines(1:52, SE_coefs$Group_1$tsa_coef, col = "red")
#lines(1:52, SE_coefs_base$Group_1$tsa_base, col = "coral", lty = 2, lwd = 1.5)
lines(1:52, SE_coefs$Group_1$aao_coef, col = "darkgreen")
#lines(1:52, SE_coefs_base$Group_1$aao_base, col = "green", lty = 2, lwd = 1.5)
segments(x0 = -1, y0 = 0, x1 = 54, y1 = 0, lty = 2)
legend("topright", inset = c(-0.3, 0),
       legend = legend1_names,
       col = legend1_col, 
       lty = legend1_lty, lwd = 2, bty = "o") #, bg = "white")
dev.off()



