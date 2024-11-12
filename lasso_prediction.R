#LASSO predictions based on the NE and SE Aus fused LASSO models

##TODO:
# try a link function for log(response), look at histograms first
# clean this code up after everything is working.
# perform CV to select lambda values (and gamma if use genlasso)

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

# grouping functions

NE_preds <- NElag_grouping(NE_laglist = NE_laglist_std, j = -19)
NE_resp <- NEresp_grouping(NEAus_mat = NEAus_mat, j = -19)

SE_preds <- SElag_grouping(SE_laglist = SE_laglist_std, j = -19)
SE_resp <- SEresp_grouping(SEAus_mat = SEAus_mat, j = -19)



# distance matrices
D <- getD1d(52)
empty_D <- matrix(rep(0, (52*51)) , ncol = 52, nrow = 51)
D1 <- cbind(D, empty_D, empty_D, empty_D)
D2 <- cbind(empty_D, D, empty_D, empty_D)
D3 <- cbind(empty_D, empty_D, D, empty_D)
D4 <- cbind(empty_D, empty_D, empty_D, D)

D_new <- rbind(D1, D2, D3, D4)

#matrix for intercept vector
D_int <- cbind(D_new, rep(0, length.out = 204)) 
D_int <- rbind( c(1, rep(0, length.out = 208)), D_int)

# NE Fusion Predictions

#TODO: create for loop and lists for the following fused lassos
## we can then adapt this to other methods (k-folds, etc) to determine our params.

NE_gamma <- 0.75

NEfuse_baselist <- list()
NEfuse_grouplist <- list()

n <- length(NE_resp)
for (i in 1:n) {
  NEresp_temp <- NE_resp[[i]]
  NEpred_temp <- as.matrix(NE_preds[[i]][ ,1:208])
}

#group 1
NE_resp1 <- NE_resp[[1]]
NE_pred1 <- as.matrix(NE_preds[[1]][ ,1:208])

NEfuse_base1 <- fusedlasso(y = NE_resp1, X = NE_pred1, D_new, gamma = 1) #with gamma = 1

NEfuse_group1 <- fusedlasso(y = NE_resp1, X = NE_pred1, D_new, gamma = NE_gamma)

#group 2 
NE_resp2 <- NE_resp[[2]]
NE_pred2 <- as.matrix(NE_preds[[2]][ ,1:208])

NEfuse_base2 <- fusedlasso(y = NE_resp2, X = NE_pred2, D_new, gamma = 1) #with gamma = 1

NEfuse_group2 <- fusedlasso(y = NE_resp2, X = NE_pred2, D_new, gamma = NE_gamma)

#group 3
NE_resp3 <- NE_resp[[3]]
NE_pred3 <- as.matrix(NE_preds[[3]][ ,1:208])

NEfuse_base3 <- fusedlasso(y = NE_resp3, X = NE_pred3, D_new, gamma = 1) #with gamma = 1

NEfuse_group3 <- fusedlasso(y = NE_resp3, X = NE_pred3, D_new, gamma = NE_gamma)

#group 4
NE_resp4 <- NE_resp[[4]]
NE_pred4 <- as.matrix(NE_preds[[4]][ ,1:208])

NEfuse_base4 <- fusedlasso(y = NE_resp4, X = NE_pred4, D_new, gamma = 1) #with gamma = 1

NEfuse_group4 <- fusedlasso(y = NE_resp4, X = NE_pred4, D_new, gamma = NE_gamma)

#group 5
NE_resp5 <- NE_resp[[5]]
NE_pred5 <- as.matrix(NE_preds[[5]][ ,1:208])

NEfuse_base5 <- fusedlasso(y = NE_resp5, X = NE_pred5, D_new, gamma = 1) #with gamma = 1

NEfuse_group5 <- fusedlasso(y = NE_resp5, X = NE_pred5, D_new, gamma = NE_gamma)

#group 6
NE_resp6 <- NE_resp[[6]]
NE_pred6 <- as.matrix(NE_preds[[6]][ ,1:208])

NEfuse_base6 <- fusedlasso(y = NE_resp6, X = NE_pred6, D_new, gamma = 1) #with gamma = 1

NEfuse_group6 <- fusedlasso(y = NE_resp6, X = NE_pred6, D_new, gamma = NE_gamma)



# Coefficient Visualizations

#TODO: create for loop to streamline this a little
## pre-define the legend vectors
## set_up lists

for(k in 1:length(NEfuse_grouplist)){
  #TODO: move below work into here
}

# group 1 coef
n <- length(NE_resp1)
p <- ncol(NE_pred1)
lambda1 <- sqrt(n * log(p))

NEgroup1_coef <- coef(NEfuse_group1, lambda = lambda1)
NEbase1_coef <- coef(NEfuse_base1, lambda = lambda1)

#extract each index
nino_coef1 <- NEgroup1_coef$beta[1:52,]
dmi_coef1 <- NEgroup1_coef$beta[53:104,]
tsa_coef1 <- NEgroup1_coef$beta[105:156,]
aao_coef1 <- NEgroup1_coef$beta[157:208,]

nino_base1 <- NEbase1_coef$beta[1:52,]
dmi_base1 <- NEbase1_coef$beta[53:104,]
tsa_base1 <- NEbase1_coef$beta[105:156,]
aao_base1 <- NEbase1_coef$beta[157:208,]

#TODO: mess with the par and output so that the legend fits well
## (check on older code for this)
par(mar = c(5, 4, 4, 8) + 0.2) #, xpd = TRUE)  # Increase the right margin (fourth value) to make space for the legend
plot(1:52, nino_coef1, type = "l", xlab = "Lag", ylab = "Coefficients",
     main = "NE Group 1", ylim = range(NEgroup1_coef$beta), 
     col = "darkmagenta")
lines(1:52, nino_base1, col = "magenta", lty = 2, lwd = 1.5)
lines(1:52, dmi_coef1, col = "blue")
lines(1:52, dmi_base1, col = "cyan", lty = 2, lwd = 1.5)
lines(1:52, tsa_coef1, col = "red")
lines(1:52, tsa_base1, col = "coral", lty = 2, lwd = 1.5)
lines(1:52, aao_coef1, col = "darkgreen")
lines(1:52, aao_base1, col = "green", lty = 2, lwd = 1.5)
abline(h = 0, lty = 2)
legend("topright", inset = c(-0.5, 0),
       legend = c("Group 1", "Group 1 - Base"),
       col = c("darkmagenta", "magenta"), 
       lty = c(1, 2), lwd = 2, bty = "o") #, bg = "white")


# Predictions


