#For model validation and visualization

#libraries
suppressMessages( library(hierNet))
suppressMessages( library(RAMP))
suppressMessages(library(glmnet)) #test ridge regression for coefs

suppressMessages( library(lubridate))
suppressMessages( library(fields)) #for some viz and other functions

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
#SE_group1_min <- hierNet(X_1, y_1, lam = cv_group$lamhat, strong = TRUE, diagonal = TRUE)

#SE Aus group 2
j <- 2
y_1 <- as.numeric(SE_resp[[j]])
X_1 <- cbind(as.matrix(SE_preds[[j]][ ,1:260]),
             as.matrix(SE_preds_q75[[j]][ ,1:104])  )
cv_group <- SEcv_new$Group_2
SE_group2_1se <- hierNet(X_1, y_1, lam = cv_group$lamhat.1se, strong = TRUE, diagonal = TRUE)
#SE_group2_min <- hierNet(X_1, y_1, lam = cv_group$lamhat, strong = TRUE, diagonal = TRUE)

#SE Aus group 3
j <- 3
y_1 <- as.numeric(SE_resp[[j]])
X_1 <- cbind(as.matrix(SE_preds[[j]][ ,1:260]),
             as.matrix(SE_preds_q75[[j]][ ,1:104])  )
cv_group <- SEcv_new2$Group_3
SE_group3_test <- hierNet(X_1, y_1, lam = cv_group$lamlist[9], strong = TRUE, diagonal = TRUE) #testing for BIC min
SE_group3_1se <- hierNet(X_1, y_1, lam = cv_group$lamhat.1se, strong = TRUE, diagonal = TRUE)
#SE_group3_min <- hierNet(X_1, y_1, lam = cv_group$lamhat, strong = TRUE, diagonal = TRUE)

#SE Aus group 4
j <- 4
y_1 <- as.numeric(SE_resp[[j]])
X_1 <- cbind(as.matrix(SE_preds[[j]][ ,1:260]),
             as.matrix(SE_preds_q75[[j]][ ,1:104])  )
cv_group <- SEcv_new2$Group_4
SE_group4_1se <- hierNet(X_1, y_1, lam = cv_group$lamhat.1se, strong = TRUE, diagonal = TRUE)
#SE_group4_min <- hierNet(X_1, y_1, lam = cv_group$lamhat, strong = TRUE, diagonal = TRUE)

#add to list
#SE_fit1_strong <- list(SE_group1_min, SE_group2_min, SE_group3_min, SE_group4_min)
SE_fit2_strong <- list(SE_group1_1se, SE_group2_1se, SE_group3_1se, SE_group4_1se)

SE_fit_test <- list(SE_group1_1se, SE_group2_1se, SE_group3_test, SE_group4_1se)


##NE Aus
#NE Aus group 1
j <- 1
y_1 <- as.numeric(NE_resp[[j]])
X_1 <- cbind(as.matrix(NE_preds[[j]][ ,1:260]),
             as.matrix(NE_preds_q75[[j]][ ,1:104])  )
cv_group <- NEcv_new$Group_1
NE_group1_1se <- hierNet(X_1, y_1, lam = cv_group$lamhat.1se, strong = TRUE, diagonal = TRUE)
#NE_group1_min <- hierNet(X_1, y_1, lam = cv_group$lamhat, strong = TRUE, diagonal = TRUE)

#NE Aus group 2
j <- 2
y_1 <- as.numeric(NE_resp[[j]])
X_1 <- cbind(as.matrix(NE_preds[[j]][ ,1:260]),
             as.matrix(NE_preds_q75[[j]][ ,1:104])  )
cv_group <- NEcv_new$Group_2
NE_group2_1se <- hierNet(X_1, y_1, lam = cv_group$lamhat.1se, strong = TRUE, diagonal = TRUE)
#NE_group2_min <- hierNet(X_1, y_1, lam = cv_group$lamhat, strong = TRUE, diagonal = TRUE)

#NE Aus group 3
j <- 3
y_1 <- as.numeric(NE_resp[[j]])
X_1 <- cbind(as.matrix(NE_preds[[j]][ ,1:260]),
             as.matrix(NE_preds_q75[[j]][ ,1:104])  )
cv_group <- NEcv_new$Group_3
NE_group3_1se <- hierNet(X_1, y_1, lam = cv_group$lamhat.1se, strong = TRUE, diagonal = TRUE)
#NE_group3_min <- hierNet(X_1, y_1, lam = cv_group$lamhat, strong = TRUE, diagonal = TRUE)

#NE Aus group 4
j <- 4
y_1 <- as.numeric(NE_resp[[j]])
X_1 <- cbind(as.matrix(NE_preds[[j]][ ,1:260]),
             as.matrix(NE_preds_q75[[j]][ ,1:104])  )
cv_group <- NEcv_new$Group_4
NE_group4_1se <- hierNet(X_1, y_1, lam = cv_group$lamhat.1se, strong = TRUE, diagonal = TRUE)
#NE_group4_min <- hierNet(X_1, y_1, lam = cv_group$lamhat, strong = TRUE, diagonal = TRUE)

#add to list
#NE_fit1_strong <- list(NE_group1_min, NE_group2_min, NE_group3_min, NE_group4_min)
NE_fit2_strong <- list(NE_group1_1se, NE_group2_1se, NE_group3_1se, NE_group4_1se)


setwd("~/CO_AUS/Aus_CO-main/Interactions")
#save(SE_fit1_strong, SE_fit2_strong, NE_fit1_strong, NE_fit2_strong, file = "hiernet_strongfits.rda")
save(SE_fit2_strong, NE_fit2_strong, file = "hiernet_strongfits.rda")

#TODO: finalize the below work (eg. refit, bic, ebic)
set.panel(2,1)
plot(NEcv_new[[3]])
plot(SEcv_new2[[1]])

#temp work:

path_group <- NEpath_new$Group_3
cv_group <- NEcv_new$Group_3

cv_group$lamlist
cv_group$lamhat.1se
cv_group$lamhat

lambda_mse <- which(cv_group$lamlist == cv_group$lamhat)
lambda_1se <- which(cv_group$lamlist == cv_group$lamhat.1se)

ebic.gamma <- 0.25 #testing gamma values in eBIC

NE3_refit1 <- list() #refit lm
NE3_refit2 <- list() #refit ridge
BIC_matrix <- matrix(NA, ncol = 5)
colnames(BIC_matrix) <- c("lambda", "lm.BIC", "ridge.BIC", "lm.eBIC", "ridge.eBIC")



for (k in 1:lambda_mse) {
  #k <- 11
  temp_bp <- path_group$bp[,k]
  temp_bn <- path_group$bn[,k]
  
  temp_th <- path_group$th[,,k]
  
  j <- 3 #group 1
  #SE Aus data
  # resp = SE_resp[[j]]
  # preds = SE_preds[[j]]
  # preds_quant = SE_preds_q75[[j]]
  
  #NE Aus data
  resp = NE_resp[[j]]
  preds = NE_preds[[j]]
  preds_quant = NE_preds_q75[[j]]  
  
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
  
  mains
  main_effect[main_terms]
  
  #length(main_terms)
  #length(interactions)
  
  #interaction effects
  interact_effect <- temp_th
  interact_terms <- which(interact_effect != 0, arr.ind = TRUE)
  
  interact_effect[interact_terms]
  
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
  #setup model for other fits
  model_string <- paste( paste(mains, collapse = " + "),
                         paste(interactions, collapse = " + "),
                         sep = " + ")
  
  model_string <- paste0("y ~ ", model_string)
  
  
  #elif
  #linear model
  data_df <- as.data.frame(cbind(y,X))
  lm_fit <- lm(formula(model_string), data = data_df)
  lm_coef <- coef(lm_fit)
  
  #ridge model
  X_df <- as.data.frame(X)
  f <- as.formula(model_string)
  X_new <- model.matrix(f, X_df)
  
  set.seed(300)
  ridge_cv <- cv.glmnet(X_new, y, alpha = 0.00, nfolds = 5)
  ridge_fit <- glmnet(X_new, y, alpha = 0.00)
  coef_pred <- predict(ridge_fit, s = ridge_cv$lambda.1se, type = "coefficients")
  par_vec <- coef_pred@i + 1 
  ridge_coef <- coef_pred[par_vec, ] #coefs from ridge
  
  lm_resid <- y - predict(lm_fit, X_df)
  ridge_resid <- y - predict(ridge_fit, X_new, s = ridge_cv$lambda.1se, type = "response")
  
  #BIC
  lm_BIC <- length(y)*log(mean(lm_resid^2)) + log(length(y))*length(lm_coef[-1])
  ridge_BIC <- length(y)*log(mean(ridge_resid^2)) + log(length(y))*length(ridge_coef[-1])
  
  #eBIC 
  #get p effective.
  p <- dim(X)[2]
  df.main <- length(main_terms)
  p_eff <- p + df.main *(df.main + 1)/2
  
  lm_eBIC <- lm_BIC + 2 * ebic.gamma * log(choose(p_eff, length(lm_coef[-1])))
  ridge_eBIC <- ridge_BIC +  2 * ebic.gamma * log(choose(p_eff, length(ridge_coef[-1])))
  
  NE3_refit1[[paste0("LamIndex_", k)]] <- lm_fit
  NE3_refit2[[paste0("LamIndex_", k)]] <- ridge_fit
  
  BIC_matrix <- rbind(BIC_matrix, c(cv_group$lamlist[k], lm_BIC, ridge_BIC, lm_eBIC, ridge_eBIC))
}

BIC_matrix <- BIC_matrix[-1, ]

lmFit_new <-NE3_refit1[[8]] 
summary(lmFit_new)

  
rm(k)

### Get coefficients for group 3

#for SE Aus Group 3
#SE_group3_1se
temp_bp <- SE_group3_1se$bp
temp_bn <- SE_group3_1se$bn

temp_th <- SE_group3_1se$th

j <- 3 #group 1
#SE Aus data
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

mains
main_effect[main_terms]

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
  
  interact_effect[interact_terms]
}



#TODO: compare these later



# predictions
NE_season <- list()
SE_season <- list()

for (k in 1:19) {
#k <- 19
  
  SEtemp_preds <- SElag_grouping(SE_laglist = SE_laglist_std, j = k)
  SEtemp_resp <- SEresp_grouping(SEAus_mat = SEAus_mat, j = k)
  SEtemp_preds_q75 <- SElag_grouping(SE_laglist = SE_laglist_q75, j = k)
  
  NEtemp_preds <- NElag_grouping(NE_laglist = NE_laglist_std, j = k)
  NEtemp_resp <- NEresp_grouping(NEAus_mat = NEAus_mat, j = k)
  NEtemp_preds_q75 <- NElag_grouping(NE_laglist = NE_laglist_q75, j = k)
  
  
  se_yhat <- list()
  se_yhat2 <- list()
  
  seresp_new <- c()
  n <- length(SEtemp_resp)
  for (i in 1:n) {
    y_temp <- as.numeric(SEtemp_resp[[i]])
    seresp_new <- c(seresp_new, y_temp)
    
    X_temp <- cbind(as.matrix(SEtemp_preds[[i]][ ,1:260]),
                 as.matrix(SEtemp_preds_q75[[i]][ ,1:104])  )
    
    #yhat <- predict(SE_fit1_strong[[i]], X_temp) #strong,lambda min
    yhat_2 <- predict(SE_fit2_strong[[i]], X_temp) #strong,lambda 1se
  
    #se_yhat[[paste0("Group_", i)]] <- yhat
    se_yhat2[[paste0("Group_", i)]] <- yhat_2
    
  }
  
  #y_sehat <- c(se_yhat[[1]], se_yhat[[2]], se_yhat[[3]], se_yhat[[4]]) 
  y_sehat2 <- c(se_yhat2[[1]], se_yhat2[[2]], se_yhat2[[3]], se_yhat2[[4]]) 
  
  
  ne_yhat <- list()
  ne_yhat2 <- list()
  
  neresp_new <- c()
  n <- length(NEtemp_resp)
  for (i in 1:n) {
    y_temp <- as.numeric(NEtemp_resp[[i]])
    neresp_new <- c(neresp_new, y_temp)
    
    X_temp <- cbind(as.matrix(NEtemp_preds[[i]][ ,1:260]),
                    as.matrix(NEtemp_preds_q75[[i]][ ,1:104])  )
    
    #yhat <- predict(NE_fit1_strong[[i]], X_temp) #strong,lambda min
    yhat_2 <- predict(NE_fit2_strong[[i]], X_temp) #strong,lambda 1se
    
    #ne_yhat[[paste0("Group_", i)]] <- yhat
    ne_yhat2[[paste0("Group_", i)]] <- yhat_2
    
  }
  
  #y_nehat <- c(ne_yhat[[1]], ne_yhat[[2]], ne_yhat[[3]], ne_yhat[[4]]) 
  y_nehat2 <- c(ne_yhat2[[1]], ne_yhat2[[2]], ne_yhat2[[3]], ne_yhat2[[4]]) 
  
  NE_season[[paste0("Season_", seasons[k])]] <- y_nehat2
  SE_season[[paste0("Season_", seasons[k])]] <- y_sehat2
}

#Figures:

#all year predictions
new_resp <- bounded_resp_df[which(bounded_resp_df$week %in% season_weeks, arr.ind = TRUE), ]
new_resp <- new_resp[-c(1:14), ]

#setup vertical year lines
unique_yr <- unique(year(new_resp$time))
year_lines <- c(unique_yr[1:length(unique_yr)])
year_lines <- paste0(year_lines, "0101")
year_lines <- as_date(year_lines)

#year_text <- year_lines
#year_text <- as_date(year_text) + months(6)
year_text <- year_lines + months(6)

xlim_val <- ymd(c("20010829", "20200401"))

temp_NE <- scale(new_resp$NE_Aus_anomaly_co, center = TRUE, scale = FALSE) 
temp_SE <- scale(new_resp$SE_Aus_anomaly_co, center = TRUE, scale = FALSE) 
temp_time <- as.Date(new_resp$time)

#prediction data:
NE_hat <- c()
SE_hat <- c()
for (i in 1:19) {
  NE_hat <- c(NE_hat, NE_season[[i]])
  SE_hat <- c(SE_hat, SE_season[[i]])
}

#region DFs
NEAUS_df <- data.frame(time = temp_time, true = temp_NE, pred = NE_hat)
SEAUS_df <- data.frame(time = temp_time, true = temp_SE, pred = SE_hat)

#test by removing every other row in the above dfs
#remove_row <- seq(1, nrow(NEAUS_df), by = 2)
#NEAUS_df <- NEAUS_df[-remove_row, ]
#SEAUS_df <- SEAUS_df[-remove_row, ]
#TODO: set up 

setwd("~/CO_AUS/Aus_CO-main/Interactions/Figures")

png(filename = "allSeasons_preds.png", width = 2200, height = 1800, res = 200)
#par(mar = c(3,5,4,1))
set.panel(2,1)
plot(x = NEAUS_df$time, y = NEAUS_df$true, pch = 18, ylim = c(min(NEAUS_df$true)-5, max(NEAUS_df$true)),
     col = "black", xlim = xlim_val, cex = 1.33, cex.main = 1.5,
     xaxt = "n",  xlab = "Year",
     ylab = "Atmospheric CO", 
     main = "NE Aus : Predictions")
points(x = NEAUS_df$time, y = NEAUS_df$pred, pch = 17, 
       col = "magenta3", cex = 1)
abline(v = year_lines[-1], lty = 2)
abline(h = 0, lty = 2)
legend("bottomleft", inset = c(0.00, 0),
       legend = c("True", "Pred"),
       col = c("black", "magenta3"),
       pch = c(18,17),
       pt.cex = c(1.33,1))
axis(side = 1, at = year_text, labels = unique_yr, tick = FALSE, cex.axis = 1)

plot(x = SEAUS_df$time, y = SEAUS_df$true, pch = 18, ylim = c(min(SEAUS_df$true)-5, max(SEAUS_df$true)),
     col = "black", xlim = xlim_val, cex = 1.33, cex.main = 1.5,
     xaxt = "n",  xlab = "Year",
     ylab = "Atmospheric CO", 
     main = "SE Aus : Predictions")
points(x = SEAUS_df$time, y = SEAUS_df$pred, pch = 17, 
       col = "magenta3", cex = 1)
abline(v = year_lines[-1], lty = 2)
abline(h = 0, lty = 2)
legend("bottomleft", inset = c(0.00, 0),
       legend = c("True", "Pred"),
       col = c("black", "magenta3"),
       pch = c(18,17),
       pt.cex = c(1.33,1))
axis(side = 1, at = year_text, labels = unique_yr, tick = FALSE, cex.axis = 1)
dev.off()


#single season (2019/2020)
setwd("~/CO_AUS/Aus_CO-main/Interactions/Figures")

x_vals <- 1:32
y_lim <- range(neresp_new, seresp_new)

png(filename = "NewSeason.png", width = 2200, height = 1800, res = 200)
set.panel(2,1)
plot(x_vals, neresp_new, type = "l", lwd = 2, ylim = y_lim, 
     ylab = "Atmospheric CO",
     xlab = "Week",  main = paste0("NE Aus: ", seasons[k]), axes = FALSE, 
     cex.main = 2.25, cex.lab = 1.5, cex.axis = 1.5)
box()
axis(2)
axis(1, at = 1:32, labels = c(season_weeks),
     las = 3, cex.axis = 1.5)
#lines(1:32, y_nehat, col = "magenta3", lwd = 2, lty = 2)
lines(1:32, y_nehat2, col = "magenta3", lwd = 2.5, lty = 1)
abline(h = 0, lty = 2)
abline(v = c(12.5, 17.5, 21.5), lty = 2, col = "red" )
legend("topright",
       legend = c("True", "Predicted"),
       col = c("black", "magenta3"),
       lty = c(1, 1),
       lwd = c(2, 2))
text(5.5, -13, "Group 1", col = "red", cex = 1)
text(15, -13, "Group 2", col = "red", cex = 1)
text(19.5, -13, "Group 3", col = "red", cex = 1)
text(24.5, -13, "Group 4", col = "red", cex = 1)

plot(x_vals, seresp_new, type = "l", lwd = 2, ylim = y_lim, 
     ylab = "Atmospheric CO",
     xlab = "Week",  main = paste0("SE Aus: ", seasons[k]), axes = FALSE, 
     cex.main = 2.25, cex.lab = 1.5, cex.axis = 1.5)
box()
axis(2)
axis(1, at = 1:32, labels = c(season_weeks),
     las = 3, cex.axis = 1.5)
#lines(1:32, y_sehat, col = "magenta3", lwd = 2, lty = 2)
lines(1:32, y_sehat2, col = "magenta3", lwd = 2.5, lty = 1)
abline(h = 0, lty = 2)
abline(v = c(7.5, 16.5, 20.5), lty = 2, col = "red" )
legend("topright",
       legend = c("True", "Predicted"),
       col = c("black", "magenta3"),
       lty = c(1, 1),
       lwd = c(2, 2))
text(3.5, -13, "Group 1", col = "red", cex = 1)
text(12, -13, "Group 2", col = "red", cex = 1)
text(18.5, -13, "Group 3", col = "red", cex = 1)
text(24.5, -13, "Group 4", col = "red", cex = 1)
dev.off()


NE_resids <- matrix(NA, ncol = 32)
SE_resids <- matrix(NA, ncol = 32)

for (k in 1:19) {
  temp_NE <- NEbase_matrix[k, ] - NE_season[[k]]
  temp_SE <- SEbase_matrix[k, ] - SE_season[[k]]
  
  NE_resids <- rbind(NE_resids, temp_NE)
  SE_resids <- rbind(SE_resids, temp_SE)
}

NE_resids <- NE_resids[-1, ]
SE_resids <- SE_resids[-1, ]

setwd("~/CO_AUS/Aus_CO-main/Interactions/Figures")

png(filename = "fullresids.png", width = 2200, height = 1800, res = 200)
set.panel(2,1)
boxplot(NE_resids, ylim = c(-20,20), ylab = "Residuals", xlab = "Week",
        main = "NE Aus : Model Residuals",axes = FALSE, pch = 20, 
        cex.main = 2.25, cex.lab = 1.5, cex.axis = 1.5)
box()
axis(2)
axis(1, at = 1:32, labels = c(season_weeks),
     las = 3, cex.axis = 1.6)
abline(h = 0, lty = 2)
abline(v = c(12.5, 17.5, 21.5), lty =2, col = "red")
text(7, -15, "Group 1", col = "red", cex =1)
text(15, -15, "Group 2", col = "red", cex = 1)
text(19.5, -15, "Group 3", col = "red", cex = 1)
text(27, -15, "Group 4", col = "red", cex = 1)

boxplot(SE_resids, ylim = c(-20,20), ylab = "Residuals", xlab = "Week",
        main = "SE Aus : Model Residuals",axes = FALSE, pch = 20, 
        cex.main = 2.25, cex.lab = 1.5, cex.axis = 1.5)
box()
axis(2)
axis(1, at = 1:32, labels = c(season_weeks),
     las = 3, cex.axis = 1.6)
abline(h = 0, lty = 2)
abline(v = c(7.5, 16.5, 20.5), lty =2, col = "red")
text(4, -15, "Group 1", col = "red", cex =1)
text(12, -15, "Group 2", col = "red", cex = 1)
text(18.5, -15, "Group 3", col = "red", cex = 1)
text(26, -15, "Group 4", col = "red", cex = 1)
dev.off()
