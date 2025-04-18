## ridge refit and model selection ##

#TODO: this will be the finalized refit coefficients,
## previous analysis work for final model selection is done elsewhere.

#libraries
suppressMessages( library(hierNet)) 
suppressMessages(library(glmnet)) 
suppressMessages( library( lubridate))
suppressMessages(library(grid))

suppressMessages( library(scales)) #for adjusting opacity

# data and functions
setwd("~/CO_AUS/Aus_CO-main/Interactions")

load( "Data/ne_data.rda")
load( "Data/se_data.rda")
load( "Data/bounded_data.rda")
load( "Data/data_matrix.rda")
load( "Data/lag_list.rda")
load( "Data/data_quantile.rda")

#functions:
source("group_functions.R") #grouping/clustering
source("refit_functions.R") #coef/refit functions

#model fits
load( "full_strong.rda") #include indicator functions
output_ind <- output

load( "main_strong.rda") #full strong models #output models
output_main <- output

load( "base_strong.rda") #strong models without 2019/2020
output_base <- output

#set up 

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

#full model data set-up
NE_preds <- NElag_grouping(NE_laglist = NE_laglist_std, j = 1:19)
NE_resp <- NEresp_grouping(NEAus_mat = NEAus_mat, j = 1:19)

SE_preds <- SElag_grouping(SE_laglist = SE_laglist_std, j = 1:19)
SE_resp <- SEresp_grouping(SEAus_mat = SEAus_mat, j = 1:19)

#quantile 75 indicator (if needed)
NE_preds_q75 <- NElag_grouping(NE_laglist = NE_laglist_q75, j = 1:19)
SE_preds_q75 <- SElag_grouping(SE_laglist = SE_laglist_q75, j = 1:19)

## refiting steps
#1. model selection via lasso AIC
#2. refit with ridge
#3. refit with ridge w/o 2019/2020 
## -compare coefficients between the same selected model

#(repeat with lasso w/o 2019/2020; variables selected with only 2000-2019)

#temp_cv <- hierNet::hierNet.cv(NEAus1$hierpath, NE_preds[[1]],  NE_resp[[1]], nfolds = 10)
#NE Aus Group 1
NEAus1 <- output_main[[1]]

NEcoefs1 <- get_coefs(NEAus1$hierpath, max_index = 50, 
                      resp = NE_resp[[1]], preds = NE_preds[[1]], preds_quant = NE_preds_q75[[1]])

plot(NEAus1$hiercv)

NEcv1 <- NEAus1$hiercv
which(NEcv1$lamlist == NEcv1$lamhat.1se)
which(NEcv1$lamlist == NEcv1$lamhat)
NEcv1$lamhat.1se
1-NEcv1$cv.err/var( NE_resp[[1]])

#length(NEcoefs1[[3]][[1]]$Main_Effect)
#length(NEcoefs1[[3]][[2]]$Interact_Effect)

NErefit1 <- refit_bic(NEAus1$hierpath, max_index = 50, ebic.gamma = 0.25, lambda.min = TRUE, AIC.c = FALSE,
                      resp = NE_resp[[1]], preds = NE_preds[[1]], preds_quant = NE_preds_q75[[1]])

NErefit1[[5]]

which.min(NErefit1[[5]][1:21, 2])
which.min(NErefit1[[5]][ , 6])

i <- 3
ridge_coef <- NErefit1[[4]][[i]][-1]
NEridge1_main <- ridge_coef
ridge_coef
length(ridge_coef)

NEcoefs1[[i]][[1]]$Main_Effect
NEcoefs1[[i]][[1]]$Coef

#NErefit1[[2]][[i]][-1] ##glm coef (fit)

#refit without 2019/2020
NE_preds <- NElag_grouping(NE_laglist = NE_laglist_std, j = -c(19))
NE_resp <- NEresp_grouping(NEAus_mat = NEAus_mat, j = -c(19))
NE_preds_q75 <- NElag_grouping(NE_laglist = NE_laglist_q75, j = -c(19))

#main model, e.g. fit var selected but change coefs

NErefit1 <- refit_bic(NEAus1$hierpath, max_index = 50, ebic.gamma = 0.25, lambda.min = TRUE, AIC.c = FALSE,
                      resp = NE_resp[[1]], preds = NE_preds[[1]], preds_quant = NE_preds_q75[[1]])

NErefit1[[5]]

#base model, e.g. model w/o 2019/2020 var select

NEAus1 <- output_base[[1]]

NEcoefs1 <- get_coefs(NEAus1$hierpath, max_index = 50, 
                      resp = NE_resp[[1]], preds = NE_preds[[1]], preds_quant = NE_preds_q75[[1]])

plot(NEAus1$hiercv)

NEcv1 <- NEAus1$hiercv
which(NEcv1$lamlist == NEcv1$lamhat.1se)
which(NEcv1$lamlist == NEcv1$lamhat)
NEcv1$lamhat.1se
1-NEcv1$cv.err/var( NE_resp[[1]])

NErefit1 <- refit_bic(NEAus1$hierpath, max_index = 50, ebic.gamma = 0.25, lambda.min = TRUE, AIC.c = FALSE,
                      resp = NE_resp[[1]], preds = NE_preds[[1]], preds_quant = NE_preds_q75[[1]])

NErefit1[[5]]

which.min(NErefit1[[5]][1:21, 2])

i <- 1
ridge_coef <- NErefit1[[4]][[i]][-1]
NEridge1_base <- ridge_coef
ridge_coef
length(ridge_coef)

#get data for all lambda values
#for() {}

#TODO: move this down after everything is setup
# Set layout
layout(matrix(c(1, 6,
                2, 6,
                3, 6,
                4, 6,
                5, 6), ncol = 2, byrow = TRUE),
       widths = c(1.5, 1), heights = c(1, 1, 1, 1, 1))

# Store info for linking
links <- list()

# --- Data Set-up --- 
#main - nino
NE1_ninolag <- c(5:12)
NE1_ninocoef <- NEridge1_main [1:8]
#base - nino (w/o 2019/2020)
NE2_ninolag <- c(7:12)
NE2_ninocoef <- NEridge1_base [1:6]
#main - tsa
NE1_tsalag <- c(14,15)
NE1_tsacoef <- NEridge1_main [9:10]


NEAus1_range <- range(NE1_ninocoef, NE2_ninocoef, NE1_tsacoef)

# --- Plot 1: Nino ---
par(mar = c(4, 4, 2, 1))
#NE Aug Group 1 

plot(NE1_ninolag, NE1_ninocoef, pch = 15, col = "green4", xlim = c(1,52), cex = 1.2,
        ylim = NEAus1_range,
        xlab = "Lag", ylab = "Coefficients")
points(NE2_ninolag, NE2_ninocoef, pch = 17, col = "cyan3", cex = 1.2)
abline(h = 0, lty = 2)

# --- Plot 3: TSA ---
par(mar = c(4, 4, 2, 1))

plot(NE1_tsalag, NE1_tsacoef, pch = 15, col = "darkorange2", xlim = c(1,52), cex = 1.2,
     ylim = NEAus1_range,
     xlab = "Lag", ylab = "Coefficients")
abline(h = 0, lty = 2)

#TODO: add in plot for coefficients


#NE Aus Group 2
NE_preds <- NElag_grouping(NE_laglist = NE_laglist_std, j = 1:19)
NE_resp <- NEresp_grouping(NEAus_mat = NEAus_mat, j = 1:19)
#quantile 75 indicator (if needed)
NE_preds_q75 <- NElag_grouping(NE_laglist = NE_laglist_q75, j = 1:19)
#main 
NEAus2 <- output_main[[2]]

NEcoefs2 <- get_coefs(NEAus2$hierpath, max_index = 50, 
                      resp = NE_resp[[2]], preds = NE_preds[[2]], preds_quant = NE_preds_q75[[2]])

plot(NEAus2$hiercv)

NEcv2 <- NEAus2$hiercv
which(NEcv2$lamlist == NEcv2$lamhat.1se)
which(NEcv2$lamlist == NEcv2$lamhat)
NEcv2$lamhat.1se
1-NEcv2$cv.err/var( NE_resp[[2]])

#length(NEcoefs2[[14]][[1]]$Main_Effect)
#length(NEcoefs2[[14]][[2]]$Interact_Effect)

NErefit2 <- refit_bic(NEAus2$hierpath, max_index = 50, ebic.gamma = 0.25, lambda.min = TRUE, AIC.c = FALSE,
                      resp = NE_resp[[2]], preds = NE_preds[[2]], preds_quant = NE_preds_q75[[2]])

NErefit2[[5]]

#plot(1:50, NErefit2[[5]][,3], type = "b", pch = 16)

AIC_index <- which.min(NErefit2[[5]][,2])
BIC_index <- which.min(NErefit2[[5]][,6])

i <- 18 #actual min at 4 (5-7)
ridge_coef <- NErefit2[[4]][[i]][-1]
lm_coef <- NErefit2[[2]][[i]][-1]
NEridge2_main <- ridge_coef
ridge_coef
lm_coef
length(ridge_coef)

lm_out <- NErefit2[[1]][[i]]


NEridge_cv <- NErefit2[[6]]
lm1se <- which( NEridge_cv[[i]]$lambda == NEridge_cv[[i]]$lambda.1se)
lmhat <- which( NEridge_cv[[i]]$lambda == NEridge_cv[[i]]$lambda.min)
NEr2 <- 1-NEridge_cv[[i]]$cvm/var( NE_resp[[2]])
NEr2[lm1se]
NEr2[lmhat]


#refit without 2019/2020
NE_preds <- NElag_grouping(NE_laglist = NE_laglist_std, j = -c(19))
NE_resp <- NEresp_grouping(NEAus_mat = NEAus_mat, j = -c(19))
NE_preds_q75 <- NElag_grouping(NE_laglist = NE_laglist_q75, j = -c(19))

#main model, e.g. fit var selected but change coefs

NEAus2 <- output_main[[2]]

NErefit2 <- refit_bic(NEAus2$hierpath, max_index = 50, ebic.gamma = 0.25, lambda.min = TRUE, AIC.c = FALSE,
                      resp = NE_resp[[2]], preds = NE_preds[[2]], preds_quant = NE_preds_q75[[2]])

NErefit2[[5]]

which.min(NErefit2[[5]][,2])
which.min(NErefit2[[5]][,6])

i <- 18 
ridge_coef <- NErefit2[[4]][[i]][-1]
NEridge2_reduce <- ridge_coef
ridge_coef
length(ridge_coef)

#base model, e.g. model w/o 2019/2020 var select

NEAus2 <- output_base[[2]]

NEcoefs2 <- get_coefs(NEAus2$hierpath, max_index = 50, 
                      resp = NE_resp[[2]], preds = NE_preds[[2]], preds_quant = NE_preds_q75[[2]])

plot(NEAus2$hiercv)

NEcv2 <- NEAus2$hiercv
which(NEcv2$lamlist == NEcv2$lamhat.1se)
which(NEcv2$lamlist == NEcv2$lamhat)
NEcv2$lamhat.1se
1-NEcv2$cv.err/var( NE_resp[[2]])

#length(NEcoefs2[[14]][[1]]$Main_Effect)
#length(NEcoefs2[[14]][[2]]$Interact_Effect)

NErefit2 <- refit_bic(NEAus2$hierpath, max_index = 50, ebic.gamma = 0.25, lambda.min = TRUE, AIC.c = FALSE,
                      resp = NE_resp[[2]], preds = NE_preds[[2]], preds_quant = NE_preds_q75[[2]])

NErefit2[[5]]

which.min(NErefit2[[5]][,2])
which.min(NErefit2[[5]][,6])

i <- 8 #actual min at 4 (5-7)
ridge_coef <- NErefit2[[4]][[i]][-1]
NEridge2_base <- ridge_coef
ridge_coef
length(ridge_coef)



#NE Aus group 2 plots

setwd("~/CO_AUS/Aus_CO-main/Interactions/Figures")

png(filename = "NE2_coeff.png", width = 3000, height = 3000, res = 300)
# Set layout
layout(matrix(c(1, 6,
                2, 6,
                3, 6,
                4, 6,
                5, 6), ncol = 2, byrow = TRUE),
       widths = c(1.5, 1), heights = c(1, 1, 1, 1, 1))

# Store info for linking
links <- list()

# --- Data Set-up --- 
#TODO: remove lag values for identical selected features (e.g. refit coefs)
#main - nino
NE1_ninolag <- c(11,15,16,25,30)
NE1_ninocoef <- NEridge2_main[1:5]
#main - nino w/o 2019/2020)
NE15_ninolag <- c(11,15,16,25,30)
NE15_ninocoef <- NEridge2_reduce[1:5]
#base - nino (w/o 2019/2020)
NE2_ninolag <- c(16,25)
NE2_ninocoef <- NEridge2_base[1:2]

#main - dmi
NE1_dmilag <- c(7,12,41,50)
NE1_dmicoef <- NEridge2_main[6:9]
#main - dmi (w/o 2019/2020)
NE15_dmilag <- c(7,12,41,50)
NE15_dmicoef <- NEridge2_reduce[6:9]
#base - dmi (w/o 2019/2020)
NE2_dmilag <- c(3)
NE2_dmicoef <- NEridge2_base[3]

#main - tsa
NE1_tsalag <- c(31,34,44,47:49)
NE1_tsacoef <- NEridge2_main[10:15]
#main - tsa (w/o 2019/2020)
NE15_tsalog <- c(31,34,44,47:49)
NE15_tsacoef <- NEridge2_reduce[10:15]
#base - tsa (w/o 2019/2020)
NE2_tsalag <- c(47, 48)
NE2_tsacoef <- NEridge2_base[4:5]

#main - sam (aao)
NE1_aaolag <- c(38, 40, 52)
NE1_aaocoef <- NEridge2_main[16:18]

#main - olr
NE1_olrlag <- c(1,2,30,32,43,52)
NE1_olrcoef <- NEridge2_main[19:24]
#main - olr (w/o 2019/2020)
NE15_olrlag <- c(1,2,30,32,43,52) 
NE15_olrcoef <- NEridge2_reduce[19:24]
#base - olr (w/o 2019/2020)
NE2_olrlag <- c(1)
NE2_olrcoef <- NEridge2_base[6]

#range of values for all coefs, for plot
NEAus2_range <- range(NE1_ninocoef, NE1_dmicoef, NE1_tsacoef, NE1_aaocoef, NE1_olrcoef,
                      NE15_ninocoef, NE15_dmicoef, NE15_tsacoef, NE15_olrcoef,
                      NE2_ninocoef, NE2_dmicoef, NE2_tsacoef, NE2_olrcoef)

# --- Plot 1: Nino ---
par(mar = c(4, 4, 2, 1))
#NE Aug Group 1 

#TODO: cahnge pch = 15 to an outlined 
plot(NE1_ninolag, NE1_ninocoef, pch = 22, 
     col = "grey4", bg =  alpha("green4",.95), cex = 1.2,
     xlim = c(1,52), 
     ylim = NEAus2_range,
     xlab = "Lag", ylab = "Coefficients")
points(NE15_ninolag, NE15_ninocoef, pch = 22, col = "black",
       bg =  alpha("chartreuse3",.65) , cex = 1.25)
points(NE2_ninolag, NE2_ninocoef, pch = 17, col = "chartreuse3", cex = 1.25)
abline(h = 0, lty = 2)
title("Nino", adj = 0)

# --- Plot 2: DMI ---
par(mar = c(4, 4, 2, 1))

plot(NE1_dmilag, NE1_dmicoef, pch = 22,
     col = "grey4",
     bg =  alpha("magenta4",.95), cex = 1.2,
     xlim = c(1,52), 
     ylim = NEAus2_range,
     xlab = "Lag", ylab = "Coefficients")
points(NE15_dmilag, NE15_dmicoef, pch = 22, col = "black",
       bg =  alpha("darkorchid2",.65) , cex = 1.25)
points(NE2_dmilag, NE2_dmicoef, pch = 17, col = "darkorchid2", cex = 1.25)
abline(h = 0, lty = 2)
title("DMI", adj = 0)

# --- Plot 3: TSA ---
par(mar = c(4, 4, 2, 1))

plot(NE1_tsalag, NE1_tsacoef, pch = 15, col = "darkorange2", xlim = c(1,52), cex = 1.2,
     ylim = NEAus2_range,
     xlab = "Lag", ylab = "Coefficients")
points(NE2_tsalag, NE2_tsacoef, pch = 17, col = "darkgoldenrod4", cex = 1.25)
abline(h = 0, lty = 2)
title("TSA", adj = 0)

# --- Plot 4: SAM ---
par(mar = c(4, 4, 2, 1))

plot(NE1_aaolag, NE1_aaocoef, pch = 15, col = "red3", xlim = c(1,52), cex = 1.2,
     ylim = NEAus2_range,
     xlab = "Lag", ylab = "Coefficients")
abline(h = 0, lty = 2)
title("SAM (AAO)", adj = 0)

# --- Plot 5: OLR ---
par(mar = c(4, 4, 2, 1))

plot(NE1_olrlag, NE1_olrcoef, pch = 15, col = "blue3", xlim = c(1,52), cex = 1.2,
     ylim = NEAus2_range,
     xlab = "Lag", ylab = "Coefficients")
points(NE2_olrlag, NE2_olrcoef, pch = 17, col = "deepskyblue3", cex = 1.25)
abline(h = 0, lty = 2)
title("OLR", adj = 0)

dev.off()


# NE Aus Group 3
#full model data set-up
NE_preds <- NElag_grouping(NE_laglist = NE_laglist_std, j = 1:19)
NE_resp <- NEresp_grouping(NEAus_mat = NEAus_mat, j = 1:19)
#quantile 75 indicator (if needed)
NE_preds_q75 <- NElag_grouping(NE_laglist = NE_laglist_q75, j = 1:19)

NEAus3 <- output_main[[3]]

NEcoefs3 <- get_coefs(NEAus3$hierpath, max_index = 50, 
                      resp = NE_resp[[3]], preds = NE_preds[[3]], preds_quant = NE_preds_q75[[3]])

plot(NEAus3$hiercv)

NEcv3 <- NEAus3$hiercv
which(NEcv3$lamlist == NEcv3$lamhat.1se)
which(NEcv3$lamlist == NEcv3$lamhat)
NEcv3$lamhat.1se

1-NEcv3$cv.err/var( NE_resp[[3]])

NErefit3 <- refit_bic(NEAus3$hierpath, max_index = 50, ebic.gamma = 0.25, lambda.min = TRUE, AIC.c = FALSE,
                      resp = NE_resp[[3]], preds = NE_preds[[3]], preds_quant = NE_preds_q75[[3]])
NErefit3[[5]]


which.min(NErefit3[[5]][1:19,2])
which.min(NErefit3[[5]][,6])

which.min(NErefit3[[5]][,9])

i <- 15
ridge_coef <- NErefit3[[4]][[i]][-1]
NEridge3_main <- ridge_coef
ridge_coef
length(ridge_coef)


#refit without 2019/2020
NE_preds <- NElag_grouping(NE_laglist = NE_laglist_std, j = -c(19))
NE_resp <- NEresp_grouping(NEAus_mat = NEAus_mat, j = -c(19))
NE_preds_q75 <- NElag_grouping(NE_laglist = NE_laglist_q75, j = -c(19))

NEAus3 <- output_base[[3]]

NEcoefs3 <- get_coefs(NEAus3$hierpath, max_index = 50, 
                      resp = NE_resp[[3]], preds = NE_preds[[3]], preds_quant = NE_preds_q75[[3]])

plot(NEAus3$hiercv)

NEcv3 <- NEAus3$hiercv
which(NEcv3$lamlist == NEcv3$lamhat.1se)
which(NEcv3$lamlist == NEcv3$lamhat)
NEcv3$lamhat.1se

1-NEcv3$cv.err/var( NE_resp[[3]])

NErefit3 <- refit_bic(NEAus3$hierpath, max_index = 50, ebic.gamma = 0.25, lambda.min = TRUE, AIC.c = FALSE,
                      resp = NE_resp[[3]], preds = NE_preds[[3]], preds_quant = NE_preds_q75[[3]])
NErefit3[[5]]

which.min(NErefit3[[5]][ ,2])
which.min(NErefit3[[5]][,6])

i <- 14
ridge_coef <- NErefit3[[4]][[i]][-1]
NEridge3_base <- ridge_coef
ridge_coef
length(ridge_coef)





# NE Aus Group 4
#full model data set-up
NE_preds <- NElag_grouping(NE_laglist = NE_laglist_std, j = 1:19)
NE_resp <- NEresp_grouping(NEAus_mat = NEAus_mat, j = 1:19)
#quantile 75 indicator (if needed)
NE_preds_q75 <- NElag_grouping(NE_laglist = NE_laglist_q75, j = 1:19)

NEAus4 <- output_main[[4]]

NEcoefs4 <- get_coefs(NEAus4$hierpath, max_index = 50, 
                      resp = NE_resp[[4]], preds = NE_preds[[4]], preds_quant = NE_preds_q75[[4]])

plot(NEAus4$hiercv)

NEcv4 <- NEAus4$hiercv
which(NEcv4$lamlist == NEcv4$lamhat.1se)
which(NEcv4$lamlist == NEcv4$lamhat)
NEcv4$lamhat.1se
1-NEcv4$cv.err/var( NE_resp[[4]])


NErefit4 <- refit_bic(NEAus4$hierpath, max_index = 50, ebic.gamma = 0.25, lambda.min = TRUE, AIC.c = FALSE,
                      resp = NE_resp[[4]], preds = NE_preds[[4]], preds_quant = NE_preds_q75[[4]])
NErefit4[[5]]

which.min(NErefit4[[5]][ ,2])
which.min(NErefit4[[5]][,6])

i <- 12
ridge_coef <- NErefit4[[4]][[i]][-1]
NEridge4_main <- ridge_coef
ridge_coef
length(ridge_coef)

#refit without 2019/2020
NE_preds <- NElag_grouping(NE_laglist = NE_laglist_std, j = -c(19))
NE_resp <- NEresp_grouping(NEAus_mat = NEAus_mat, j = -c(19))
NE_preds_q75 <- NElag_grouping(NE_laglist = NE_laglist_q75, j = -c(19))

NEAus4 <- output_base[[4]]

NEcoefs4 <- get_coefs(NEAus4$hierpath, max_index = 50, 
                      resp = NE_resp[[4]], preds = NE_preds[[4]], preds_quant = NE_preds_q75[[4]])

plot(NEAus4$hiercv)

NEcv4 <- NEAus4$hiercv
which(NEcv4$lamlist == NEcv4$lamhat.1se)
which(NEcv4$lamlist == NEcv4$lamhat)
NEcv4$lamhat.1se
1-NEcv4$cv.err/var( NE_resp[[4]])


NErefit4 <- refit_bic(NEAus4$hierpath, max_index = 50, ebic.gamma = 0.25, lambda.min = TRUE, AIC.c = FALSE,
                      resp = NE_resp[[4]], preds = NE_preds[[4]], preds_quant = NE_preds_q75[[4]])
NErefit4[[5]]

which.min(NErefit4[[5]][ ,2])
which.min(NErefit4[[5]][,6])

i <- 12
ridge_coef <- NErefit4[[4]][[i]][-1]
NEridge4_base <- ridge_coef
ridge_coef
length(ridge_coef)


## SE Aus Group 1

#full model data set-up
SE_preds <- SElag_grouping(SE_laglist = SE_laglist_std, j = 1:19)
SE_resp <- SEresp_grouping(SEAus_mat = SEAus_mat, j = 1:19)
SE_preds_q75 <- SElag_grouping(SE_laglist = SE_laglist_q75, j = 1:19)

SEAus1 <- output_main[[5]]

SEcoefs1 <- get_coefs(SEAus1$hierpath, max_index = 50, 
                      resp = SE_resp[[1]], preds = SE_preds[[1]], preds_quant = SE_preds_q75[[1]])

plot(SEAus1$hiercv)

SEcv1 <- SEAus1$hiercv
which(SEcv1$lamlist == SEcv1$lamhat.1se)
which(SEcv1$lamlist == SEcv1$lamhat)
SEcv1$lamhat.1se
1-SEcv1$cv.err/var( SE_resp[[1]])


SErefit1 <- refit_bic(SEAus1$hierpath, max_index = 50, ebic.gamma = 0.25, lambda.min = FALSE, AIC.c = FALSE,
                      resp = SE_resp[[1]], preds = SE_preds[[1]], preds_quant = SE_preds_q75[[1]])

SErefit1[[5]]

which.min(SErefit1[[5]][,2])
which.min(SErefit1[[5]][,6])

i <- 7
ridge_coef <- SErefit1[[4]][[i]][-1]
SEridge1_main <- ridge_coef
ridge_coef
length(ridge_coef)

##refit without 2019/2020
SE_preds <- SElag_grouping(SE_laglist = SE_laglist_std, j = -c(19))
SE_resp <- SEresp_grouping(SEAus_mat = SEAus_mat, j = -c(19))
SE_preds_q75 <- SElag_grouping(SE_laglist = SE_laglist_q75, j = -c(19))

SEAus1 <- output_base[[5]]

SEcoefs1 <- get_coefs(SEAus1$hierpath, max_index = 50, 
                      resp = SE_resp[[1]], preds = SE_preds[[1]], preds_quant = SE_preds_q75[[1]])

plot(SEAus1$hiercv)

SEcv1 <- SEAus1$hiercv
which(SEcv1$lamlist == SEcv1$lamhat.1se)
which(SEcv1$lamlist == SEcv1$lamhat)
SEcv1$lamhat.1se
1-SEcv1$cv.err/var( SE_resp[[1]])

SErefit1 <- refit_bic(SEAus1$hierpath, max_index = 50, ebic.gamma = 0.25, lambda.min = FALSE, AIC.c = FALSE,
                      resp = SE_resp[[1]], preds = SE_preds[[1]], preds_quant = SE_preds_q75[[1]])

SErefit1[[5]]

which.min(SErefit1[[5]][,2])
which.min(SErefit1[[5]][,6])

i <- 7
ridge_coef <- SErefit1[[4]][[i]][-1]
SEridge1_base <- ridge_coef
ridge_coef
length(ridge_coef)



## SE Aus Group 2
SE_preds <- SElag_grouping(SE_laglist = SE_laglist_std, j = 1:19)
SE_resp <- SEresp_grouping(SEAus_mat = SEAus_mat, j = 1:19)
SE_preds_q75 <- SElag_grouping(SE_laglist = SE_laglist_q75, j = 1:19)

SEAus2 <- output_main[[6]]

SEcoefs2 <- get_coefs(SEAus2$hierpath, max_index = 50, 
                      resp = SE_resp[[2]], preds = SE_preds[[2]], preds_quant = SE_preds_q75[[2]])

plot(SEAus2$hiercv)

SEcv2 <- SEAus2$hiercv
which(SEcv2$lamlist == SEcv2$lamhat.1se)
which(SEcv2$lamlist == SEcv2$lamhat)
SEcv2$lamhat.1se
1-SEcv2$cv.err/var( SE_resp[[2]])


SErefit2 <- refit_bic(SEAus2$hierpath, max_index = 50, ebic.gamma = 0.25, lambda.min = FALSE, AIC.c = FALSE,
                      resp = SE_resp[[2]], preds = SE_preds[[2]], preds_quant = SE_preds_q75[[2]])

SErefit2[[5]]

which.min(SErefit2[[5]][,2])
which.min(SErefit2[[5]][,6])

i <- 7
ridge_coef <- SErefit2[[4]][[i]][-1]
SEridge2_main <- ridge_coef
ridge_coef
length(ridge_coef)

##refit without 2019/2020
SE_preds <- SElag_grouping(SE_laglist = SE_laglist_std, j = -c(19))
SE_resp <- SEresp_grouping(SEAus_mat = SEAus_mat, j = -c(19))
SE_preds_q75 <- SElag_grouping(SE_laglist = SE_laglist_q75, j = -c(19))

SEAus2 <- output_base[[6]]

SEcoefs2 <- get_coefs(SEAus2$hierpath, max_index = 50, 
                      resp = SE_resp[[2]], preds = SE_preds[[2]], preds_quant = SE_preds_q75[[2]])

plot(SEAus2$hiercv)

SEcv2 <- SEAus2$hiercv
which(SEcv2$lamlist == SEcv2$lamhat.1se)
which(SEcv2$lamlist == SEcv2$lamhat)
SEcv2$lamhat.1se
1-SEcv2$cv.err/var( SE_resp[[2]])


SErefit2 <- refit_bic(SEAus2$hierpath, max_index = 50, ebic.gamma = 0.25, lambda.min = FALSE, AIC.c = FALSE,
                      resp = SE_resp[[2]], preds = SE_preds[[2]], preds_quant = SE_preds_q75[[2]])

SErefit2[[5]]

which.min(SErefit2[[5]][,2])
which.min(SErefit2[[5]][,6])

i <- 8
ridge_coef <- SErefit2[[4]][[i]][-1]
SEridge2_base <- ridge_coef
ridge_coef
length(ridge_coef)



## SE Aus Group 3
SE_preds <- SElag_grouping(SE_laglist = SE_laglist_std, j = 1:19)
SE_resp <- SEresp_grouping(SEAus_mat = SEAus_mat, j = 1:19)
SE_preds_q75 <- SElag_grouping(SE_laglist = SE_laglist_q75, j = 1:19)

SEAus3 <- output_main[[7]]

SEcoefs3 <- get_coefs(SEAus3$hierpath, max_index = 50, 
                      resp = SE_resp[[3]], preds = SE_preds[[3]], preds_quant = SE_preds_q75[[3]])

plot(SEAus3$hiercv)

SEcv3 <- SEAus3$hiercv
which(SEcv3$lamlist == SEcv3$lamhat.1se)
which(SEcv3$lamlist == SEcv3$lamhat)
SEcv3$lamhat.1se
1-SEcv3$cv.err/var( SE_resp[[3]])

SErefit3 <- refit_bic(SEAus3$hierpath, max_index = 50, ebic.gamma = 0.50, lambda.min = TRUE, AIC.c = FALSE,
                      resp = SE_resp[[3]], preds = SE_preds[[3]], preds_quant = SE_preds_q75[[3]])
SErefit3[[5]]

which.min(SErefit3[[5]][,2])
which.min(SErefit3[[5]][,6])

i <- 17
ridge_coef <- SErefit3[[4]][[i]][-1]
SEridge3_main <- ridge_coef
ridge_coef
length(ridge_coef)

##refit without 2019/2020
SE_preds <- SElag_grouping(SE_laglist = SE_laglist_std, j = -c(19))
SE_resp <- SEresp_grouping(SEAus_mat = SEAus_mat, j = -c(19))
SE_preds_q75 <- SElag_grouping(SE_laglist = SE_laglist_q75, j = -c(19))


SEAus3 <- output_base[[7]]

SEcoefs3 <- get_coefs(SEAus3$hierpath, max_index = 50, 
                      resp = SE_resp[[3]], preds = SE_preds[[3]], preds_quant = SE_preds_q75[[3]])

plot(SEAus3$hiercv)

SEcv3 <- SEAus3$hiercv
which(SEcv3$lamlist == SEcv3$lamhat.1se)
which(SEcv3$lamlist == SEcv3$lamhat)
SEcv3$lamhat.1se
1-SEcv3$cv.err/var( SE_resp[[3]])

SErefit3 <- refit_bic(SEAus3$hierpath, max_index = 50, ebic.gamma = 0.50, lambda.min = TRUE, AIC.c = FALSE,
                      resp = SE_resp[[3]], preds = SE_preds[[3]], preds_quant = SE_preds_q75[[3]])
SErefit3[[5]]

which.min(SErefit3[[5]][,2])
which.min(SErefit3[[5]][,6])

i <- 4
ridge_coef <- SErefit3[[4]][[i]][-1]
SEridge3_base <- ridge_coef
ridge_coef
length(ridge_coef)


## SE Aus Group 4
SE_preds <- SElag_grouping(SE_laglist = SE_laglist_std, j = 1:19)
SE_resp <- SEresp_grouping(SEAus_mat = SEAus_mat, j = 1:19)
SE_preds_q75 <- SElag_grouping(SE_laglist = SE_laglist_q75, j = 1:19)

SEAus4 <- output_main[[8]]

SEcoefs4 <- get_coefs(SEAus4$hierpath, max_index = 50, 
                      resp = SE_resp[[4]], preds = SE_preds[[4]], preds_quant = SE_preds_q75[[4]])

plot(SEAus4$hiercv)

SEcv4 <- SEAus4$hiercv
which(SEcv4$lamlist == SEcv4$lamhat.1se)
which(SEcv4$lamlist == SEcv4$lamhat)
SEcv4$lamhat.1se
#1-SEcv4$cv.err/var( SE_resp[[4]])


SErefit4 <- refit_bic(SEAus4$hierpath, max_index = 50, ebic.gamma = 0.25, lambda.min = TRUE, AIC.c = FALSE,
                      resp = SE_resp[[4]], preds =  SE_preds[[4]], preds_quant = SE_preds_q75[[4]])

SErefit4[[5]]

which.min(SErefit4[[5]][,2])
which.min(SErefit4[[5]][,6])

i <- 18
ridge_coef <- SErefit4[[4]][[i]][-1]
SEridge4_main <- ridge_coef
ridge_coef
length(ridge_coef)

##refit without 2019/2020
SE_preds <- SElag_grouping(SE_laglist = SE_laglist_std, j = -c(19))
SE_resp <- SEresp_grouping(SEAus_mat = SEAus_mat, j = -c(19))
SE_preds_q75 <- SElag_grouping(SE_laglist = SE_laglist_q75, j = -c(19))


SEAus4 <- output_base[[8]]

SEcoefs4 <- get_coefs(SEAus4$hierpath, max_index = 50, 
                      resp = SE_resp[[4]], preds = SE_preds[[4]], preds_quant = SE_preds_q75[[4]])

plot(SEAus4$hiercv)

SEcv4 <- SEAus4$hiercv
which(SEcv4$lamlist == SEcv4$lamhat.1se)
which(SEcv4$lamlist == SEcv4$lamhat)
SEcv4$lamhat.1se
#1-SEcv4$cv.err/var( SE_resp[[4]])


SErefit4 <- refit_bic(SEAus4$hierpath, max_index = 50, ebic.gamma = 0.25, lambda.min = TRUE, AIC.c = FALSE,
                      resp = SE_resp[[4]], preds =  SE_preds[[4]], preds_quant = SE_preds_q75[[4]])

SErefit4[[5]]

which.min(SErefit4[[5]][,2])
which.min(SErefit4[[5]][,6])

i <- 4
ridge_coef <- SErefit4[[4]][[i]][-1]
SEridge4_base <- ridge_coef
ridge_coef
length(ridge_coef)


