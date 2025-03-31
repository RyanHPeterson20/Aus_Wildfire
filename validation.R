#For model validation and visualization

#libraries
suppressMessages( library(hierNet))
suppressMessages( library(RAMP))
suppressMessages(library(glmnet)) #test ridge regression for coefs
suppressMessages( library(lubridate))

#data 
setwd("~/CO_AUS/Aus_CO-main/Interactions")

load( "Data/ne_data.rda")
load( "Data/se_data.rda")
load( "Data/bounded_data.rda")
load( "Data/data_matrix.rda")
load( "Data/lag_list.rda")
load( "Data/data_quantile.rda")

#functions
source("group_functions.R") #grouping/clustering

#models
load( "hiernet1_group.rda") #weak/4group

load( "SEAus_hiernet_temp2.rda") #strong/se/4group
load( "SEAus_hiernet_temp.rda") #strong/se/4group
load( "NEAus_hiernet_strong.rda") #strong/ne/4group

#TODO: load 
#load( "hiernet_group_split.rda") #strong/2group


## set-up
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

#group response matrices
#NEAus
NEAus_1 <- NEbase_matrix[ ,1:12] #early season
NEAus_2 <- NEbase_matrix[ ,13:17] #primary NE fire season
NEAus_3 <- NEbase_matrix[ ,18:21] #feedback from SE
NEAus_4 <- NEbase_matrix[ ,22:32] #late season

#SEAus
SEAus_1 <- SEbase_matrix[ ,1:7] #early season
SEAus_2 <- SEbase_matrix[ ,8:16] #feedback from NE
SEAus_3 <- SEbase_matrix[ ,17:20] #primary SE fire season
SEAus_4 <- SEbase_matrix[ ,21:32] #late season

NEAus_mat <- list(NEAus_1, NEAus_2, NEAus_3, NEAus_4)
SEAus_mat <- list(SEAus_1, SEAus_2, SEAus_3, SEAus_4)

rm(NEAus_1, NEAus_2, NEAus_3, NEAus_4, SEAus_1, SEAus_2, SEAus_3, SEAus_4)


#data setup (4 group)
#full data sets (training)
NE_preds <- NElag_grouping(NE_laglist = NE_laglist_std, j = 1:19)
NE_resp <- NEresp_grouping(NEAus_mat = NEAus_mat, j = 1:19)

SE_preds <- SElag_grouping(SE_laglist = SE_laglist_std, j = 1:19)
SE_resp <- SEresp_grouping(SEAus_mat = SEAus_mat, j = 1:19)

#quantile 75 indicator
NE_preds_q75 <- NElag_grouping(NE_laglist = NE_laglist_q75, j = 1:19)
SE_preds_q75 <- SElag_grouping(SE_laglist = SE_laglist_q75, j = 1:19)

#strong hier fits

##SE Aus
#SE Aus group 1
j <- 1
y_1 <- as.numeric(SE_resp[[j]])
X_1 <- cbind(as.matrix(SE_preds[[j]][ ,1:260]),
             as.matrix(SE_preds_q75[[j]][ ,1:104])  )
cv_group <- SEcv_new$Group_1
SE_group1_1se <- hierNet(X_1, y_1, lam = cv_group$lamhat.1se, strong = TRUE, diagonal = TRUE)
SE_group1_min <- hierNet(X_1, y_1, lam = cv_group$lamhat, strong = TRUE, diagonal = TRUE)

#SE Aus group 2
j <- 2
y_1 <- as.numeric(SE_resp[[j]])
X_1 <- cbind(as.matrix(SE_preds[[j]][ ,1:260]),
             as.matrix(SE_preds_q75[[j]][ ,1:104])  )
cv_group <- SEcv_new$Group_2
SE_group2_1se <- hierNet(X_1, y_1, lam = cv_group$lamhat.1se, strong = TRUE, diagonal = TRUE)
SE_group2_min <- hierNet(X_1, y_1, lam = cv_group$lamhat, strong = TRUE, diagonal = TRUE)

#SE Aus group 3
j <- 3
y_1 <- as.numeric(SE_resp[[j]])
X_1 <- cbind(as.matrix(SE_preds[[j]][ ,1:260]),
             as.matrix(SE_preds_q75[[j]][ ,1:104])  )
cv_group <- SEcv_new2$Group_3
SE_group3_test <- hierNet(X_1, y_1, lam = cv_group$lamlist[9], strong = TRUE, diagonal = TRUE) #testing for BIC min
SE_group3_1se <- hierNet(X_1, y_1, lam = cv_group$lamhat.1se, strong = TRUE, diagonal = TRUE)
SE_group3_min <- hierNet(X_1, y_1, lam = cv_group$lamhat, strong = TRUE, diagonal = TRUE)

#SE Aus group 4
j <- 4
y_1 <- as.numeric(SE_resp[[j]])
X_1 <- cbind(as.matrix(SE_preds[[j]][ ,1:260]),
             as.matrix(SE_preds_q75[[j]][ ,1:104])  )
cv_group <- SEcv_new2$Group_4
SE_group4_1se <- hierNet(X_1, y_1, lam = cv_group$lamhat.1se, strong = TRUE, diagonal = TRUE)
SE_group4_min <- hierNet(X_1, y_1, lam = cv_group$lamhat, strong = TRUE, diagonal = TRUE)

#add to list
SE_fit1_strong <- list(SE_group1_min, SE_group2_min, SE_group3_min, SE_group4_min)
SE_fit2_strong <- list(SE_group1_1se, SE_group2_1se, SE_group3_1se, SE_group4_1se)

SE_fit_test <- list(SE_group1_1se, SE_group2_1se, SE_group3_test, SE_group4_1se)

#TODO: finish up NE Aus


#TODO: finalize the below work (eg. refit, bic, ebic)
plot(SEcv_new2[[1]])

SEcv_new2$Group_3$lamlist

path_group <- SEpath_new2$Group_3
cv_group <- SEcv_new2$Group_3

cv_group$lamhat.1se
cv_group$lamhat

lambda_mse <- which(cv_group$lamlist == cv_group$lamhat)
lambda_1se <- which(cv_group$lamlist == cv_group$lamhat.1se)

#for (k in 1:lambda_mse) {
  k <- 9
  cv_group$lamlist[k]
  
  temp_bp <- path_group$bp[,k]
  temp_bn <- path_group$bn[,k]
  
  temp_th <- path_group$th[,,k]
  
  j <- 1 #group 1
  resp = SE_resp[[j]]
  preds = SE_preds[[j]]
  preds_quant = SE_preds_q75[[j]]
  
  
  #internals
  
  y = as.numeric(resp)
  X = cbind(as.matrix(preds[ ,1:260]),
            as.matrix(preds_quant[ ,1:104])  )
  
  
  coef_names <- c(colnames(preds), paste0( "IND_", colnames(preds_quant[ ,1:104]) ) )
  colnames(X) <- coef_names
  
  #TODO: add in method for looping through different lambda
  
  #main effects
  main_effect <- temp_bp - temp_bn
  main_terms <- which(main_effect != 0, arr.ind = TRUE)
  mains <- colnames(X)[main_terms]
  
  #length(main_terms)
  #length(interactions)
  
  #interaction effects
  interact_effect <- temp_th
  interact_terms <- which(interact_effect != 0, arr.ind = TRUE)
  
  #TODO: add elif for interactions
  if (nrow(interact_terms) != 0) {
    interact_names <- matrix(coef_names[interact_terms], nrow = nrow(interact_terms), ncol = ncol(interact_terms))
    
    interactions <-  paste0(interact_names[,2], ":", interact_names[,1])
    for (i in 1:nrow(interact_names)) {
      interact_subterms <- interact_names[i, ]
      if (interact_subterms[1] != interact_subterms[2]) {
        #condition for interactions
        this_term <- paste(interact_subterms, collapse = ":")
      } else {
        #condition for squared terms
        this_term <- paste0("I(", interact_subterms[1], "^2)")
      }
      
      interactions[i] <- this_term
    }
    interactions
  }
  mains
#}

rm(k)

#TODO: compare these later



# predictions
k <- 19

SEtemp_preds <- SElag_grouping(SE_laglist = SE_laglist_std, j = k)

SEtemp_resp <- SEresp_grouping(SEAus_mat = SEAus_mat, j = k)

SEtemp_preds_q75 <- SElag_grouping(SE_laglist = SE_laglist_q75, j = k)


se_yhat <- list()
se_yhat2 <- list()
seresp_new <- c()
n <- length(SEtemp_resp)
for (i in 1:n) {
  y_temp <- as.numeric(SEtemp_resp[[i]])
  seresp_new <- c(seresp_new, y_temp)
  
  
  X_temp <- cbind(as.matrix(SEtemp_preds[[i]][ ,1:260]),
               as.matrix(SEtemp_preds_q75[[i]][ ,1:104])  )
  
  yhat <- predict(SE_fit1_strong[[i]], X_temp) #strong,lambda min
  yhat_2 <- predict(SE_fit2_strong[[i]], X_temp) #strong,lambda 1se
  
  yhat_3 <- predict(SE_fit_test[[i]], X_temp) #weak,lambda 1se
  
  se_yhat[[paste0("Group_", i)]] <- yhat
  se_yhat2[[paste0("Group_", i)]] <- yhat_2
  se_yhat3[[paste0("Group_", i)]] <- yhat_3
}

y_sehat <- c(se_yhat[[1]], se_yhat[[2]], se_yhat[[3]], se_yhat[[4]]) 
y_sehat2 <- c(se_yhat2[[1]], se_yhat2[[2]], se_yhat2[[3]], se_yhat2[[4]]) 
y_sehat3 <- c(se_yhat3[[1]], se_yhat3[[2]], se_yhat3[[3]], se_yhat3[[4]]) 

x_vals <- 1:32

#temp plot
#TODO: update with other features (eg. from AGU viz (and others))
plot(x_vals, seresp_new, type = "l", lwd = 2, ylim = range(seresp_new, y_sehat), 
     ylab = "Atmospheric CO",
     xlab = "",  main = paste0("SE Aus: ", seasons[k]), axes = FALSE, 
     cex.main = 2.25, cex.lab = 1.5, cex.axis = 1.5)
box()
axis(2)
axis(1, at = 1:32, labels = c(paste0("Week ", season_weeks)),
     las = 3, cex.axis = 1.6)
lines(1:32, y_sehat2, col = "magenta3", lwd = 2, lty = 2)
lines(1:32, y_sehat3, col = "magenta2", lwd = 2, lty = 1)

