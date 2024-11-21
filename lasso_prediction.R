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

#TODO:
## we can then adapt this to other methods (k-folds, etc) to determine our params.

#TODO: check to see if I can pass these fusedlasso objects to a predict() function

#NE_gamma <- 0.67 #original gamma
NE_gamma <- 0.85

NEfuse_baselist <- list()
NEfuse_grouplist <- list()

n <- length(NE_resp)
for (i in 1:n) {
  NEresp_temp <- NE_resp[[i]]
  NEpred_temp <- as.matrix(NE_preds[[i]][ ,1:208])
  
  NEbase_temp <- fusedlasso(y = NEresp_temp, X = NEpred_temp, D_new, gamma = 1) #with gamma = 1
  NEgroup_temp <- fusedlasso(y = NEresp_temp, X = NEpred_temp, D_new, gamma = NE_gamma)
  
  NEfuse_baselist[[paste0("Group_", i)]] <- NEbase_temp
  NEfuse_grouplist[[paste0("Group_", i)]] <- NEgroup_temp
}

# SE Fusion Predictions

#SE_gamma <- 0.67 #original gamma
SE_gamma <- 0.85


SEfuse_baselist <- list()
SEfuse_grouplist <- list()

n <- length(SE_resp)
for (i in 1:n) {
  SEresp_temp <- SE_resp[[i]]
  SEpred_temp <- as.matrix(SE_preds[[i]][ ,1:208])
  
  SEbase_temp <- fusedlasso(y = SEresp_temp, X = SEpred_temp, D_new, gamma = 1) #with gamma = 1
  SEgroup_temp <- fusedlasso(y = SEresp_temp, X = SEpred_temp, D_new, gamma = SE_gamma)
  
  SEfuse_baselist[[paste0("Group_", i)]] <- SEbase_temp
  SEfuse_grouplist[[paste0("Group_", i)]] <- SEgroup_temp
}


# Coefficient Visualizations

## pre-define the legend vectors
legend_names <- c("Nino", "Nino - Base", "DMI", "DMI - Base", 
                  "TSA", "TSA - Base", "AAO", "AAO - Base")
legend_col <- c("darkmagenta", "magenta", "blue", "cyan",
                "red", "coral", "darkgreen", "green")
legend_lty <- rep(c(1,2), length.out = 8)


legend1_names <- c("Nino",  "DMI", "TSA", "AAO")
legend1_col <- c("darkmagenta", "blue", "red", "darkgreen")
legend1_lty <- rep(1, length.out = 4)

#NE Visualizations
NE_coefs <- list()
NE_coefs_base <- list()

#full coef list
NEfuse_coefs <- list()

NE_range <- c()
NEbase_range <- c()

for(k in 1:length(NEfuse_grouplist)){

  #TODO: get better defs of lambda 
  n <- length(NE_resp[[k]])
  p <- ncol(NE_preds[[k]])
  lambda1 <- sqrt(n * log(p))
  
  NEgroup_coef <- coef(NEfuse_grouplist[[k]], lambda = lambda1)
  NEbase_coef <- coef(NEfuse_baselist[[k]], lambda = lambda1)
  
  NEfuse <- cbind(NEgroup_coef$beta, NEbase_coef$beta)
  colnames(NEfuse) <- c("NEgroup", "NEbase")
  NEfuse_coefs[[paste0("Group_", k)]] <- as.data.frame(NEfuse)
  
  NE_range <- c(NE_range, range(NEgroup_coef$beta))
  NEbase_range <- c(NEbase_range, range(NEbase_coef$beta))
  
  #extract each index
  nino_coef <- NEgroup_coef$beta[1:52,]
  dmi_coef <- NEgroup_coef$beta[53:104,]
  tsa_coef <- NEgroup_coef$beta[105:156,]
  aao_coef <- NEgroup_coef$beta[157:208,]
  
  NE_coefs[[paste0("Group_", k)]] <- data.frame(nino_coef, dmi_coef, tsa_coef, aao_coef)
  
  nino_base <- NEbase_coef$beta[1:52,]
  dmi_base <- NEbase_coef$beta[53:104,]
  tsa_base <- NEbase_coef$beta[105:156,]
  aao_base <- NEbase_coef$beta[157:208,]
  
  NE_coefs_base[[paste0("Group_", k)]] <- data.frame(nino_base, dmi_base, tsa_base, aao_base)
}



#NE group 1
png(filename = "NE_group1.png", width = 2400, height = 1800, res = 300)
par(mar = c(5, 4, 4, 9) + 0.2, xpd = TRUE)  # Increase the right margin (fourth value) to make space for the legend
plot(1:52, NE_coefs$Group_1$nino_coef, type = "l", xlab = "Lag", ylab = "Coefficients",
     main = "NE Group 1", ylim = range(NE_range, NEbase_range), 
     col = "darkmagenta")
#lines(1:52, NE_coefs_base$Group_1$nino_base, col = "magenta", lty = 2, lwd = 1.5)
lines(1:52, NE_coefs$Group_1$dmi_coef, col = "blue")
#lines(1:52, NE_coefs_base$Group_1$dmi_base, col = "cyan", lty = 2, lwd = 1.5)
lines(1:52, NE_coefs$Group_1$tsa_coef, col = "red")
#lines(1:52, NE_coefs_base$Group_1$tsa_base, col = "coral", lty = 2, lwd = 1.5)
lines(1:52, NE_coefs$Group_1$aao_coef, col = "darkgreen")
#lines(1:52, NE_coefs_base$Group_1$aao_base, col = "green", lty = 2, lwd = 1.5)
segments(x0 = -1, y0 = 0, x1 = 54, y1 = 0, lty = 2)
legend("topright", inset = c(-0.3, 0),
       legend = legend1_names,
       col = legend1_col, 
       lty = legend1_lty, lwd = 2, bty = "o") #, bg = "white")
dev.off()


#NE group 2
png(filename = "NE_group2.png", width = 2400, height = 1800, res = 300)
par(mar = c(5, 4, 4, 9) + 0.2, xpd = TRUE)  # Increase the right margin (fourth value) to make space for the legend
plot(1:52, NE_coefs$Group_2$nino_coef, type = "l", xlab = "Lag", ylab = "Coefficients",
     main = "NE Group 2", ylim = range(NE_range, NEbase_range), 
     col = "darkmagenta")
#lines(1:52, NE_coefs_base$Group_2$nino_base, col = "magenta", lty = 2, lwd = 1.5)
lines(1:52, NE_coefs$Group_2$dmi_coef, col = "blue")
#lines(1:52, NE_coefs_base$Group_2$dmi_base, col = "cyan", lty = 2, lwd = 1.5)
lines(1:52, NE_coefs$Group_2$tsa_coef, col = "red")
#lines(1:52, NE_coefs_base$Group_2$tsa_base, col = "coral", lty = 2, lwd = 1.5)
lines(1:52, NE_coefs$Group_2$aao_coef, col = "darkgreen")
#lines(1:52, NE_coefs_base$Group_2$aao_base, col = "green", lty = 2, lwd = 1.5)
segments(x0 = -1, y0 = 0, x1 = 54, y1 = 0, lty = 2)
legend("topright", inset = c(-0.3, 0),
       legend = legend1_names,
       col = legend1_col, 
       lty = legend1_lty, lwd = 2, bty = "o") #, bg = "white")
dev.off()


#NE group 3
png(filename = "NE_group3.png", width = 2400, height = 1800, res = 300)
par(mar = c(5, 4, 4, 9) + 0.2, xpd = TRUE)  # Increase the right margin (fourth value) to make space for the legend
plot(1:52, NE_coefs$Group_3$nino_coef, type = "l", xlab = "Lag", ylab = "Coefficients",
     main = "NE Group 3", ylim = range(NE_range, NEbase_range), 
     col = "darkmagenta")
#lines(1:52, NE_coefs_base$Group_3$nino_base, col = "magenta", lty = 2, lwd = 1.5)
lines(1:52, NE_coefs$Group_3$dmi_coef, col = "blue")
#lines(1:52, NE_coefs_base$Group_3$dmi_base, col = "cyan", lty = 2, lwd = 1.5)
lines(1:52, NE_coefs$Group_3$tsa_coef, col = "red")
#lines(1:52, NE_coefs_base$Group_3$tsa_base, col = "coral", lty = 2, lwd = 1.5)
lines(1:52, NE_coefs$Group_3$aao_coef, col = "darkgreen")
#lines(1:52, NE_coefs_base$Group_3$aao_base, col = "green", lty = 2, lwd = 1.5)
segments(x0 = -1, y0 = 0, x1 = 54, y1 = 0, lty = 2)
legend("topright", inset = c(-0.3, 0),
       legend = legend1_names,
       col = legend1_col, 
       lty = legend1_lty, lwd = 2, bty = "o") #, bg = "white")
dev.off()

#NE group 4
png(filename = "NE_group4.png", width = 2400, height = 1800, res = 300)
par(mar = c(5, 4, 4, 9) + 0.2, xpd = TRUE)  # Increase the right margin (fourth value) to make space for the legend
plot(1:52, NE_coefs$Group_4$nino_coef, type = "l", xlab = "Lag", ylab = "Coefficients",
     main = "NE Group 4", ylim = range(NE_range, NEbase_range), 
     col = "darkmagenta")
#lines(1:52, NE_coefs_base$Group_4$nino_base, col = "magenta", lty = 2, lwd = 1.5)
lines(1:52, NE_coefs$Group_4$dmi_coef, col = "blue")
#lines(1:52, NE_coefs_base$Group_4$dmi_base, col = "cyan", lty = 2, lwd = 1.5)
lines(1:52, NE_coefs$Group_4$tsa_coef, col = "red")
#lines(1:52, NE_coefs_base$Group_4$tsa_base, col = "coral", lty = 2, lwd = 1.5)
lines(1:52, NE_coefs$Group_4$aao_coef, col = "darkgreen")
#lines(1:52, NE_coefs_base$Group_4$aao_base, col = "green", lty = 2, lwd = 1.5)
segments(x0 = -1, y0 = 0, x1 = 54, y1 = 0, lty = 2)
legend("topright", inset = c(-0.3, 0),
       legend = legend1_names,
       col = legend1_col, 
       lty = legend1_lty, lwd = 2, bty = "o") #, bg = "white")
dev.off()


#NE group 5
png(filename = "NE_group5.png", width = 2400, height = 1800, res = 300)
par(mar = c(5, 4, 4, 9) + 0.2, xpd = TRUE)  # Increase the right margin (fourth value) to make space for the legend
plot(1:52, NE_coefs$Group_5$nino_coef, type = "l", xlab = "Lag", ylab = "Coefficients",
     main = "NE Group 5", ylim = range(NE_range, NEbase_range), 
     col = "darkmagenta")
#lines(1:52, NE_coefs_base$Group_5$nino_base, col = "magenta", lty = 2, lwd = 1.5)
lines(1:52, NE_coefs$Group_5$dmi_coef, col = "blue")
#lines(1:52, NE_coefs_base$Group_5$dmi_base, col = "cyan", lty = 2, lwd = 1.5)
lines(1:52, NE_coefs$Group_5$tsa_coef, col = "red")
#lines(1:52, NE_coefs_base$Group_5$tsa_base, col = "coral", lty = 2, lwd = 1.5)
lines(1:52, NE_coefs$Group_5$aao_coef, col = "darkgreen")
#lines(1:52, NE_coefs_base$Group_5$aao_base, col = "green", lty = 2, lwd = 1.5)
segments(x0 = -1, y0 = 0, x1 = 54, y1 = 0, lty = 2)
legend("topright", inset = c(-0.3, 0),
       legend = legend1_names,
       col = legend1_col, 
       lty = legend1_lty, lwd = 2, bty = "o") #, bg = "white")
dev.off()


#NE group 6
png(filename = "NE_group6.png", width = 2400, height = 1800, res = 300)
par(mar = c(5, 4, 4, 9) + 0.2, xpd = TRUE)  # Increase the right margin (fourth value) to make space for the legend
plot(1:52, NE_coefs$Group_6$nino_coef, type = "l", xlab = "Lag", ylab = "Coefficients",
     main = "NE Group 6", ylim = range(NE_range, NEbase_range), 
     col = "darkmagenta")
#lines(1:52, NE_coefs_base$Group_6$nino_base, col = "magenta", lty = 2, lwd = 1.5)
lines(1:52, NE_coefs$Group_6$dmi_coef, col = "blue")
#lines(1:52, NE_coefs_base$Group_6$dmi_base, col = "cyan", lty = 2, lwd = 1.5)
lines(1:52, NE_coefs$Group_6$tsa_coef, col = "red")
#lines(1:52, NE_coefs_base$Group_6$tsa_base, col = "coral", lty = 2, lwd = 1.5)
lines(1:52, NE_coefs$Group_6$aao_coef, col = "darkgreen")
#lines(1:52, NE_coefs_base$Group_6$aao_base, col = "green", lty = 2, lwd = 1.5)
segments(x0 = -1, y0 = 0, x1 = 54, y1 = 0, lty = 2)
legend("topright", inset = c(-0.3, 0),
       legend = legend1_names,
       col = legend1_col, 
       lty = legend1_lty, lwd = 2, bty = "o") #, bg = "white")
dev.off()


#SE Visualizations
SE_coefs <- list()
SE_coefs_base <- list()

#full coef list
SEfuse_coefs <- list()

SE_range <- c()
SEbase_range <- c()

for(k in 1:length(SEfuse_grouplist)){
  #TODO: move below work into here
  n <- length(SE_resp[[k]])
  p <- ncol(SE_preds[[k]])
  lambda1 <- sqrt(n * log(p))
  
  SEgroup_coef <- coef(SEfuse_grouplist[[k]], lambda = lambda1)
  SEbase_coef <- coef(SEfuse_baselist[[k]], lambda = lambda1)
  
  SEfuse <- cbind(SEgroup_coef$beta, SEbase_coef$beta)
  colnames(SEfuse) <- c("SEgroup", "SEbase")
  SEfuse_coefs[[paste0("Group_", k)]] <- as.data.frame(SEfuse)
  
  SE_range <- c(SE_range, range(SEgroup_coef$beta))
  SEbase_range <- c(SEbase_range, range(SEbase_coef$beta))
  
  #extract each index
  nino_coef <- SEgroup_coef$beta[1:52,]
  dmi_coef <- SEgroup_coef$beta[53:104,]
  tsa_coef <- SEgroup_coef$beta[105:156,]
  aao_coef <- SEgroup_coef$beta[157:208,]
  
  SE_coefs[[paste0("Group_", k)]] <- data.frame(nino_coef, dmi_coef, tsa_coef, aao_coef)
  
  nino_base <- SEbase_coef$beta[1:52,]
  dmi_base <- SEbase_coef$beta[53:104,]
  tsa_base <- SEbase_coef$beta[105:156,]
  aao_base <- SEbase_coef$beta[157:208,]
  
  SE_coefs_base[[paste0("Group_", k)]] <- data.frame(nino_base, dmi_base, tsa_base, aao_base)
}



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


#SE group 2
png(filename = "SE_group2.png", width = 2400, height = 1800, res = 300)
par(mar = c(5, 4, 4, 9) + 0.2, xpd = TRUE)  # Increase the right margin (fourth value) to make space for the legend
plot(1:52, SE_coefs$Group_2$nino_coef, type = "l", xlab = "Lag", ylab = "Coefficients",
     main = "SE Group 2", ylim = range(SE_range, SEbase_range), 
     col = "darkmagenta")
#lines(1:52, SE_coefs_base$Group_2$nino_base, col = "magenta", lty = 2, lwd = 1.5)
lines(1:52, SE_coefs$Group_2$dmi_coef, col = "blue")
#lines(1:52, SE_coefs_base$Group_2$dmi_base, col = "cyan", lty = 2, lwd = 1.5)
lines(1:52, SE_coefs$Group_2$tsa_coef, col = "red")
#lines(1:52, SE_coefs_base$Group_2$tsa_base, col = "coral", lty = 2, lwd = 1.5)
lines(1:52, SE_coefs$Group_2$aao_coef, col = "darkgreen")
#lines(1:52, SE_coefs_base$Group_2$aao_base, col = "green", lty = 2, lwd = 1.5)
segments(x0 = -1, y0 = 0, x1 = 54, y1 = 0, lty = 2)
legend("topright", inset = c(-0.3, 0),
       legend = legend1_names,
       col = legend1_col, 
       lty = legend1_lty, lwd = 2, bty = "o") #, bg = "white")
dev.off()



#SE group 3
png(filename = "SE_group3.png", width = 2400, height = 1800, res = 300)
par(mar = c(5, 4, 4, 9) + 0.2, xpd = TRUE)  # Increase the right margin (fourth value) to make space for the legend
plot(1:52, SE_coefs$Group_3$nino_coef, type = "l", xlab = "Lag", ylab = "Coefficients",
     main = "SE Group 3", ylim = range(SE_range, SEbase_range), 
     col = "darkmagenta")
#lines(1:52, SE_coefs_base$Group_3$nino_base, col = "magenta", lty = 2, lwd = 1.5)
lines(1:52, SE_coefs$Group_3$dmi_coef, col = "blue")
#lines(1:52, SE_coefs_base$Group_3$dmi_base, col = "cyan", lty = 2, lwd = 1.5)
lines(1:52, SE_coefs$Group_3$tsa_coef, col = "red")
#lines(1:52, SE_coefs_base$Group_3$tsa_base, col = "coral", lty = 2, lwd = 1.5)
lines(1:52, SE_coefs$Group_3$aao_coef, col = "darkgreen")
#lines(1:52, SE_coefs_base$Group_3$aao_base, col = "green", lty = 2, lwd = 1.5)
segments(x0 = -1, y0 = 0, x1 = 54, y1 = 0, lty = 2)
legend("topright", inset = c(-0.3, 0),
       legend = legend1_names,
       col = legend1_col, 
       lty = legend1_lty, lwd = 2, bty = "o") #, bg = "white")
dev.off()


#SE group 4
png(filename = "SE_group4.png", width = 2400, height = 1800, res = 300)
par(mar = c(5, 4, 4, 9) + 0.2, xpd = TRUE)  # Increase the right margin (fourth value) to make space for the legend
plot(1:52, SE_coefs$Group_4$nino_coef, type = "l", xlab = "Lag", ylab = "Coefficients",
     main = "SE Group 4", ylim = range(SE_range, SEbase_range), 
     col = "darkmagenta")
#lines(1:52, SE_coefs_base$Group_4$nino_base, col = "magenta", lty = 2, lwd = 1.5)
lines(1:52, SE_coefs$Group_4$dmi_coef, col = "blue")
#lines(1:52, SE_coefs_base$Group_4$dmi_base, col = "cyan", lty = 2, lwd = 1.5)
lines(1:52, SE_coefs$Group_4$tsa_coef, col = "red")
#lines(1:52, SE_coefs_base$Group_4$tsa_base, col = "coral", lty = 2, lwd = 1.5)
lines(1:52, SE_coefs$Group_4$aao_coef, col = "darkgreen")
#lines(1:52, SE_coefs_base$Group_4$aao_base, col = "green", lty = 2, lwd = 1.5)
segments(x0 = -1, y0 = 0, x1 = 54, y1 = 0, lty = 2)
legend("topright", inset = c(-0.3, 0),
       legend = legend1_names,
       col = legend1_col, 
       lty = legend1_lty, lwd = 2, bty = "o") #, bg = "white")
dev.off()


#SE group 5
png(filename = "SE_group5.png", width = 2400, height = 1800, res = 300)
par(mar = c(5, 4, 4, 9) + 0.2, xpd = TRUE)  # Increase the right margin (fourth value) to make space for the legend
plot(1:52, SE_coefs$Group_5$nino_coef, type = "l", xlab = "Lag", ylab = "Coefficients",
     main = "SE Group 5", ylim = range(SE_range, SEbase_range), 
     col = "darkmagenta")
#lines(1:52, SE_coefs_base$Group_5$nino_base, col = "magenta", lty = 2, lwd = 1.5)
lines(1:52, SE_coefs$Group_5$dmi_coef, col = "blue")
#lines(1:52, SE_coefs_base$Group_5$dmi_base, col = "cyan", lty = 2, lwd = 1.5)
lines(1:52, SE_coefs$Group_5$tsa_coef, col = "red")
#lines(1:52, SE_coefs_base$Group_5$tsa_base, col = "coral", lty = 2, lwd = 1.5)
lines(1:52, SE_coefs$Group_5$aao_coef, col = "darkgreen")
#lines(1:52, SE_coefs_base$Group_5$aao_base, col = "green", lty = 2, lwd = 1.5)
segments(x0 = -1, y0 = 0, x1 = 54, y1 = 0, lty = 2)
#abline(h = 0, lty = 2)
legend("topright", inset = c(-0.3, 0),
       legend = legend1_names,
       col = legend1_col, 
       lty = legend1_lty, lwd = 2, bty = "o") #, bg = "white")
dev.off()


# Predictions

#TODO: make the following sections of predictions more flexible with for loops
#TODO: add in weekly grouping, that way we can check if a week is in a response group

#NE groups (move up later)
NEgroup1 <- season_weeks[1:3]
NEgroup2 <- season_weeks[4:8]
NEgroup3 <- season_weeks[9:12]
NEgroup4 <- season_weeks[13:17]
NEgroup5 <- season_weeks[18:21]
NEgroup6 <- season_weeks[22:32]

#SE groups
SEgroup1 <- season_weeks[1:3]
SEgroup2 <- season_weeks[4:7]
SEgroup3 <- season_weeks[8:16]
SEgroup4 <- season_weeks[17:20]
SEgroup5 <- season_weeks[21:32]

#NE in-sample predictions

#NE group 1
NEgroup1_list <- list()

ne_resp1 <- NEAus_mat[[1]]
ne_beta1 <- NEfuse_coefs$Group_1$NEgroup
nebase_beta1 <- NEfuse_coefs$Group_1$NEbase

for (k in 1:3) {
  ne_week <- NE_laglist_std[[k]]
  
  week_preds <- as.data.frame(matrix(NA, ncol = 5))
  colnames(week_preds) <- c("True", "Pred", "Resid", "Pred_Base", "Resid_Base")
  
  for (i in 1:18) {
    temp_preds <- as.matrix(ne_week[i, 3:210])
    temp_resp <- ne_resp1[i, k]
    
    pred_temp <- temp_preds %*% ne_beta1
    temp_resid <- temp_resp - pred_temp
    
    pred_base_temp <- temp_preds %*% nebase_beta1
    temp_resid_base <- temp_resp - pred_base_temp
    
    week_preds <- rbind(week_preds, c(temp_resp, pred_temp, temp_resid, pred_base_temp, temp_resid_base))
  }
  
  week_preds <- week_preds[-1, ]
  
  NEgroup1_list[[paste0("Week_", season_weeks[k])]] <- week_preds
}


#NE group 2
NEgroup2_list <- list()

ne_resp2 <- NEAus_mat[[2]]
ne_beta2 <- NEfuse_coefs$Group_2$NEgroup
nebase_beta2 <- NEfuse_coefs$Group_2$NEbase

for (k in 4:8) {
  ne_week <- NE_laglist_std[[k]]
  
  week_preds <- as.data.frame(matrix(NA, ncol = 5))
  colnames(week_preds) <- c("True", "Pred", "Resid", "Pred_Base", "Resid_Base")
  
  week_k <- which(NEgroup2 == season_weeks[k])
  for (i in 1:18) {
    temp_preds <- as.matrix(ne_week[i, 3:210])
    temp_resp <- ne_resp2[i, week_k]
    
    pred_temp <- temp_preds %*% ne_beta2
    temp_resid <- temp_resp - pred_temp
    
    pred_base_temp <- temp_preds %*% nebase_beta2
    temp_resid_base <- temp_resp - pred_base_temp
    
    week_preds <- rbind(week_preds, c(temp_resp, pred_temp, temp_resid, pred_base_temp, temp_resid_base))
  }
  
  week_preds <- week_preds[-1, ]
  
  NEgroup2_list[[paste0("Week_", season_weeks[k])]] <- week_preds
}


#NE group 3
NEgroup3_list <- list()

ne_resp3 <- NEAus_mat[[3]]
ne_beta3 <- NEfuse_coefs$Group_3$NEgroup
nebase_beta3 <- NEfuse_coefs$Group_3$NEbase

for (k in 9:12) {
  ne_week <- NE_laglist_std[[k]]
  
  week_preds <- as.data.frame(matrix(NA, ncol = 5))
  colnames(week_preds) <- c("True", "Pred", "Resid", "Pred_Base", "Resid_Base")
  
  week_k <- which(NEgroup3 == season_weeks[k])
  for (i in 1:18) {
    temp_preds <- as.matrix(ne_week[i, 3:210])
    temp_resp <- ne_resp3[i, week_k]
    
    pred_temp <- temp_preds %*% ne_beta3
    temp_resid <- temp_resp - pred_temp
    
    pred_base_temp <- temp_preds %*% nebase_beta3
    temp_resid_base <- temp_resp - pred_base_temp
    
    week_preds <- rbind(week_preds, c(temp_resp, pred_temp, temp_resid, pred_base_temp, temp_resid_base))
  }
  
  week_preds <- week_preds[-1, ]
  
  NEgroup3_list[[paste0("Week_", season_weeks[k])]] <- week_preds
}


#NE group 4
NEgroup4_list <- list()

ne_resp4 <- NEAus_mat[[4]]
ne_beta4 <- NEfuse_coefs$Group_4$NEgroup
nebase_beta4 <- NEfuse_coefs$Group_4$NEbase

for (k in 13:17) {
  ne_week <- NE_laglist_std[[k]]
  
  week_preds <- as.data.frame(matrix(NA, ncol = 5))
  colnames(week_preds) <- c("True", "Pred", "Resid", "Pred_Base", "Resid_Base")

  week_k <- which(NEgroup4 == season_weeks[k])  
  for (i in 1:18) {
    temp_preds <- as.matrix(ne_week[i, 3:210])
    temp_resp <- ne_resp4[i, week_k]
    
    pred_temp <- temp_preds %*% ne_beta4
    temp_resid <- temp_resp - pred_temp
    
    pred_base_temp <- temp_preds %*% nebase_beta4
    temp_resid_base <- temp_resp - pred_base_temp
    
    week_preds <- rbind(week_preds, c(temp_resp, pred_temp, temp_resid, pred_base_temp, temp_resid_base))
  }
  
  week_preds <- week_preds[-1, ]
  
  NEgroup4_list[[paste0("Week_", season_weeks[k])]] <- week_preds
}


#NE group 5
NEgroup5_list <- list()

ne_resp5 <- NEAus_mat[[5]]
ne_beta5 <- NEfuse_coefs$Group_5$NEgroup
nebase_beta5 <- NEfuse_coefs$Group_5$NEbase

for (k in 18:21) {
  ne_week <- NE_laglist_std[[k]]
  
  week_preds <- as.data.frame(matrix(NA, ncol = 5))
  colnames(week_preds) <- c("True", "Pred", "Resid", "Pred_Base", "Resid_Base")
  
  week_k <- which(NEgroup5 == season_weeks[k]) 
  for (i in 1:18) {
    temp_preds <- as.matrix(ne_week[i, 3:210])
    temp_resp <- ne_resp5[i, week_k]
    
    pred_temp <- temp_preds %*% ne_beta5
    temp_resid <- temp_resp - pred_temp
    
    pred_base_temp <- temp_preds %*% nebase_beta5
    temp_resid_base <- temp_resp - pred_base_temp
    
    week_preds <- rbind(week_preds, c(temp_resp, pred_temp, temp_resid, pred_base_temp, temp_resid_base))
  }
  
  week_preds <- week_preds[-1, ]
  
  NEgroup5_list[[paste0("Week_", season_weeks[k])]] <- week_preds
}


#NE group 6
NEgroup6_list <- list()

ne_resp6 <- NEAus_mat[[6]]
ne_beta6 <- NEfuse_coefs$Group_6$NEgroup
nebase_beta6 <- NEfuse_coefs$Group_6$NEbase

for (k in 22:32) {
  ne_week <- NE_laglist_std[[k]]
  
  week_preds <- as.data.frame(matrix(NA, ncol = 5))
  colnames(week_preds) <- c("True", "Pred", "Resid", "Pred_Base", "Resid_Base")
  
  week_k <- which(NEgroup6 == season_weeks[k]) 
  for (i in 1:18) {
    temp_preds <- as.matrix(ne_week[i, 3:210])
    temp_resp <- ne_resp6[i, week_k]
    
    pred_temp <- temp_preds %*% ne_beta6
    temp_resid <- temp_resp - pred_temp
    
    pred_base_temp <- temp_preds %*% nebase_beta6
    temp_resid_base <- temp_resp - pred_base_temp
    
    week_preds <- rbind(week_preds, c(temp_resp, pred_temp, temp_resid, pred_base_temp, temp_resid_base))
  }
  
  week_preds <- week_preds[-1, ]
  
  NEgroup6_list[[paste0("Week_", season_weeks[k])]] <- week_preds
}

#SE in-sample Predictions


#SE group 1
SEgroup1_list <- list()

se_resp1 <- SEAus_mat[[1]]
se_beta1 <- SEfuse_coefs$Group_1$SEgroup
sebase_beta1 <- SEfuse_coefs$Group_1$SEbase

for (k in 1:3) {
  se_week <- SE_laglist_std[[k]]
  
  week_preds <- as.data.frame(matrix(NA, ncol = 5))
  colnames(week_preds) <- c("True", "Pred", "Resid", "Pred_Base", "Resid_Base")
  
  for (i in 1:18) {
    temp_preds <- as.matrix(se_week[i, 3:210])
    temp_resp <- se_resp1[i, k]
    
    pred_temp <- temp_preds %*% se_beta1
    temp_resid <- temp_resp - pred_temp
    
    pred_base_temp <- temp_preds %*% sebase_beta1
    temp_resid_base <- temp_resp - pred_base_temp
    
    week_preds <- rbind(week_preds, c(temp_resp, pred_temp, temp_resid, pred_base_temp, temp_resid_base))
  }
  
  week_preds <- week_preds[-1, ]
  
  SEgroup1_list[[paste0("Week_", season_weeks[k])]] <- week_preds
}


#SE group 2
SEgroup2_list <- list()

se_resp2 <- SEAus_mat[[2]]
se_beta2 <- SEfuse_coefs$Group_2$SEgroup
sebase_beta2 <- SEfuse_coefs$Group_2$SEbase

for (k in 4:7) {
  se_week <- SE_laglist_std[[k]]
  
  week_preds <- as.data.frame(matrix(NA, ncol = 5))
  colnames(week_preds) <- c("True", "Pred", "Resid", "Pred_Base", "Resid_Base")
  
  week_k <- which(SEgroup2 == season_weeks[k])
  for (i in 1:18) {
    temp_preds <- as.matrix(se_week[i, 3:210])
    temp_resp <- se_resp2[i, week_k]
    
    pred_temp <- temp_preds %*% se_beta2
    temp_resid <- temp_resp - pred_temp
    
    pred_base_temp <- temp_preds %*% sebase_beta2
    temp_resid_base <- temp_resp - pred_base_temp
    
    week_preds <- rbind(week_preds, c(temp_resp, pred_temp, temp_resid, pred_base_temp, temp_resid_base))
  }
  
  week_preds <- week_preds[-1, ]
  
  SEgroup2_list[[paste0("Week_", season_weeks[k])]] <- week_preds
}


#SE group 3
SEgroup3_list <- list()

se_resp3 <- SEAus_mat[[3]]
se_beta3 <- SEfuse_coefs$Group_3$SEgroup
sebase_beta3 <- SEfuse_coefs$Group_3$SEbase


for (k in 8:16) {
  se_week <- SE_laglist_std[[k]]
  
  week_preds <- as.data.frame(matrix(NA, ncol = 5))
  colnames(week_preds) <- c("True", "Pred", "Resid", "Pred_Base", "Resid_Base")
  
  week_k <- which(SEgroup3 == season_weeks[k])
  for (i in 1:18) {
    temp_preds <- as.matrix(se_week[i, 3:210])
    temp_resp <- se_resp3[i, week_k]
    
    pred_temp <- temp_preds %*% se_beta3
    temp_resid <- temp_resp - pred_temp
    
    pred_base_temp <- temp_preds %*% sebase_beta3
    temp_resid_base <- temp_resp - pred_base_temp
    
    week_preds <- rbind(week_preds, c(temp_resp, pred_temp, temp_resid, pred_base_temp, temp_resid_base))
  }
  
  week_preds <- week_preds[-1, ]
  
  SEgroup3_list[[paste0("Week_", season_weeks[k])]] <- week_preds
}


#SE group 4
SEgroup4_list <- list()

se_resp4 <- SEAus_mat[[4]]
se_beta4 <- SEfuse_coefs$Group_4$SEgroup
sebase_beta4 <- SEfuse_coefs$Group_4$SEbase

for (k in 17:20) {
  se_week <- SE_laglist_std[[k]]
  
  week_preds <- as.data.frame(matrix(NA, ncol = 5))
  colnames(week_preds) <- c("True", "Pred", "Resid", "Pred_Base", "Resid_Base")
  
  week_k <- which(SEgroup4 == season_weeks[k])
  for (i in 1:18) {
    temp_preds <- as.matrix(se_week[i, 3:210])
    temp_resp <- se_resp4[i, week_k]
    
    pred_temp <- temp_preds %*% se_beta4
    temp_resid <- temp_resp - pred_temp
    
    pred_base_temp <- temp_preds %*% sebase_beta4
    temp_resid_base <- temp_resp - pred_base_temp
    
    week_preds <- rbind(week_preds, c(temp_resp, pred_temp, temp_resid, pred_base_temp, temp_resid_base))
  }
  
  week_preds <- week_preds[-1, ]
  
  SEgroup4_list[[paste0("Week_", season_weeks[k])]] <- week_preds
}


#SE group 5
SEgroup5_list <- list()

se_resp5 <- SEAus_mat[[5]]
se_beta5 <- SEfuse_coefs$Group_5$SEgroup
sebase_beta5 <- SEfuse_coefs$Group_5$SEbase


for (k in 21:32) {
  se_week <- SE_laglist_std[[k]]
  
  week_preds <- as.data.frame(matrix(NA, ncol = 5))
  colnames(week_preds) <- c("True", "Pred", "Resid", "Pred_Base", "Resid_Base")
  
  week_k <- which(SEgroup5 == season_weeks[k])
  for (i in 1:18) {
    temp_preds <- as.matrix(se_week[i, 3:210])
    temp_resp <- se_resp5[i, week_k]
    
    pred_temp <- temp_preds %*% se_beta5
    temp_resid <- temp_resp - pred_temp
    
    pred_base_temp <- temp_preds %*% sebase_beta5
    temp_resid_base <- temp_resp - pred_base_temp
    
    week_preds <- rbind(week_preds, c(temp_resp, pred_temp, temp_resid, pred_base_temp, temp_resid_base))
  }
  
  week_preds <- week_preds[-1, ]
  
  SEgroup5_list[[paste0("Week_", season_weeks[k])]] <- week_preds
}



NEgroup_list <- c(NEgroup1_list, NEgroup2_list, NEgroup3_list, NEgroup4_list, NEgroup5_list, NEgroup6_list)
SEgroup_list <- c(SEgroup1_list, SEgroup2_list, SEgroup3_list, SEgroup4_list, SEgroup5_list)

#TODO: export this to mess with elsewhere
#setwd("~/CO_AUS/Aus_CO-main")

#save(NEgroup_list, SEgroup_list, NEfuse_coefs, SEfuse_coefs, file = "base_preds.rda")

#residual histograms

#TODO: add in SE residuals when they are done

residNE1_box_df <- data.frame(values = NA, group = NA)
residSE1_box_df <- data.frame(values = NA, group = NA)

for (i in 1:16) {
  temp_vals <- NEgroup_list[[i]]$Resid
  temp_group <- rep(season_weeks[i], length(NEgroup_list[[i]]$Resid))
  residNE1_box_df <- rbind(residNE1_box_df, list(temp_vals, temp_group))
  
  temp_vals1 <- SEgroup_list[[i]]$Resid
  temp_group1 <- rep(season_weeks[i], length(SEgroup_list[[i]]$Resid))
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
  temp_vals <- NEgroup_list[[i]]$Resid
  temp_group <- rep(season_weeks[i], length(NEgroup_list[[i]]$Resid))
  residNE2_box_df <- rbind(residNE2_box_df, list(temp_vals, temp_group))
  
  temp_vals1 <- SEgroup_list[[i]]$Resid
  temp_group1 <- rep(season_weeks[i], length(SEgroup_list[[i]]$Resid))
  residSE2_box_df <- rbind(residSE2_box_df, list(temp_vals1, temp_group1))
}

residNE2_box_df <- residNE2_box_df[-1, ]
residNE2_box_df$group <- as.factor(residNE2_box_df$group)
residNE2_box_df$values <- as.numeric(residNE2_box_df$values)

residSE2_box_df <- residSE2_box_df[-1, ]
residSE2_box_df$group <- as.factor(residSE2_box_df$group)
residSE2_box_df$values <- as.numeric(residSE2_box_df$values)

resid_lim <- range(residNE1_box_df$values, residNE2_box_df$values)

#resid_min <- min(residNE1_box_df$values, residNE2_box_df$values, 
#                 residSE1_box_df$values, residSE2_box_df$values)

#resid_max <- max(residNE1_box_df$values, residNE2_box_df$values, 
#                 residSE1_box_df$values, residSE2_box_df$values)

#resid_lim <- c(resid_min, resid_max)

residNEfull_box_df <- rbind(residNE1_box_df, residNE2_box_df)
residSEfull_box_df <- rbind(residSE1_box_df, residSE2_box_df)


setwd("~/CO_AUS/Aus_CO-main/Figures_Lasso")

#TODO:update using resid_lim from the simplified linear model
resid_lim <- c(-20, 20)

png("NEresidsd_newgamma.png", width = 800, height = 600, res = 100)
par(mar = c(8, 4, 4, 2) + 0.1)
boxplot(values ~ group, data = residNEfull_box_df, ylim = resid_lim, ylab = "Residuals", xlab = "",
        main = "NE Aus Lasso Model Residuals",axes = FALSE, pch = 20)
box()
axis(2)
axis(1, at = 1:32, labels = c(paste0("Week ", season_weeks)),
     las = 3)
abline(h = 0, lty = 2)
abline(v = c(3.5, 8.5, 12.5, 17.5, 21.5), lty =2, col = "red")
text(1.5, -15, "Group 1", col = "red", cex = 0.75)
text(6, -15, "Group 2", col = "red", cex = 0.75)
text(10.5, -15, "Group 3", col = "red", cex = 0.75)
text(15, -15, "Group 4", col = "red", cex = 0.75)
text(19.5, -15, "Group 5", col = "red", cex = 0.75)
text(25, -15, "Group 6", col = "red", cex = 0.75)
dev.off()


png("SEresidsd_newgamma.png", width = 800, height = 600, res = 100)
par(mar = c(8, 4, 4, 2) + 0.1)
boxplot(values ~ group, data = residSEfull_box_df, ylim = resid_lim, ylab = "Residuals", xlab = "",
        main = "SE Aus Lasso Model Residuals", axes = FALSE, pch = 20)
box()
axis(2)
axis(1, at = 1:32, labels = c(paste0("Week ", season_weeks)),
     las = 3)
abline(h = 0, lty = 2)
abline(v = c(3.5, 7.5, 16.5, 20.5), lty =2, col = "red")
text(1.5, -15, "Group 1", col = "red", cex = 0.75)
text(5.5, -15, "Group 2", col = "red", cex = 0.75)
text(12, -15, "Group 3", col = "red", cex = 0.75)
text(18.5, -15, "Group 4", col = "red", cex = 0.75)
text(24, -15, "Group 5", col = "red", cex = 0.75)
dev.off()



