
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
source("predictionTest_functions.R") #predtest functions 
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

## 
